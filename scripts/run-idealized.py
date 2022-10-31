import argparse 
import datetime as dt
import numpy as np

DAYS = 24 * 60 * 60

def roots_of_unity_upper(n):
    """(2n)-th roots of unity (positive imaginary only)"""
    return np.exp(1j * np.pi * (np.arange(n) + 0.5) / n)


class ETDRK4Diagonal:
    """Exponential time differencing 4-th oder Runge-Kutta scheme

    Supports diagonal linear operators only.

    Based on the implementation of the Cox and Matthews (2002) scheme by Kassam
    and Trefethen (2005) using contour integrals to evaluate the phi-functions.
    """

    def __init__(self, step, linear_diag, nonlinear, n_contour=16):
        self.linear_diag = np.asarray(linear_diag)
        self.nonlinear = nonlinear
        # Number of quadrature nodes for contour integral
        self._n_contour = n_contour
        self.step = step
        self._prepare()

    def _prepare(self):
        h = self.step
        L = self.linear_diag
        A = h*L
        self._E  = np.exp(    A)
        self._E2 = np.exp(0.5*A)
        # Contour integral quadrature nodes (along complex circle)
        r = roots_of_unity_upper(self._n_contour)
        # Translate nodes so circles are centered around the points of evaluation
        LR = A[:,None] + r[None,:]
        # Evaluate contour integrals
        # For a, b, c coefficients: Q = (e^Lh/2 - I) / L
        self._Q = h * np.real(np.mean((np.exp(LR/2) - 1) / LR, axis=1))
        # Alpha, beta, gamma coefficients
        self._α = h * np.real(np.mean((-4 -   LR         + np.exp(LR)*( 4 - 3*LR + (LR**2))) / (LR**3), axis=1))
        self._β = h * np.real(np.mean(( 2 +   LR         + np.exp(LR)*(-2 +   LR          )) / (LR**3), axis=1))
        self._γ = h * np.real(np.mean((-4 - 3*LR - LR**2 + np.exp(LR)*( 4 -   LR          )) / (LR**3), axis=1))

    def advance(self, u, t):
        N  = self.nonlinear
        Q  = self._Q
        E  = self._E
        E2 = self._E2
        # Evaluate the nonlinear term
        Nu = N(u, t)
        # Coefficient a and associated nonlinear term
        a  = E2*u + Q*Nu
        Na = N(a, t + self.step*0.5)
        # Coefficient b and associated nonlinear term
        b  = E2*u + Q*Na
        Nb = N(b, t + self.step*0.5)
        # Coefficient c and associated nonlinear term
        c  = E2*a + Q*(2*Nb - Nu)
        Nc = N(c, t + self.step)
        # Combine
        return E*u + Nu*self._α + 2*(Na+Nb)*self._β + Nc*self._γ

    def run(self, u, nsteps, t, keep_every=0):
        n = 0 # step counter
        us, ns = [u], [n]
        while n <= nsteps:
            u = self.advance(u, t)
            t += self.step
            n += 1
            if keep_every > 0 and n % keep_every == 0:
                us.append(u)
                ns.append(n)
        # Zero-valued keep_every means no intermediate steps are kept, only
        # initial conditions and end state
        if keep_every <= 0:
            us.append(u)
            ns.append(n)
        # Create time array from counted steps and return data numpyified
        ts = np.asarray(ns) * self.step
        ts += t - n*self.step
        us = np.asarray(us)
        return ts, us


def hann(x, x0, width):
    x = np.asarray(x)
    y = (x - x0) / width * np.pi
    out = np.zeros(y.shape)
    sup = (-0.5*np.pi <= y) & (y <= 0.5*np.pi)
    out[sup] = np.cos(y[sup])**2
    return out


parser = argparse.ArgumentParser()
parser.add_argument("outfile", type=str)
parser.add_argument("--nens", type=int, default=25)
parser.add_argument("--alpha", type=float, default=0.4, help="wave-mean flow interaction parameter")
parser.add_argument("--urefcg", type=float, default=60., help="background state group velocity")
parser.add_argument("--forcing-peak", type=float, default=-2.)

if __name__ == "__main__":
    args = parser.parse_args()

    import xarray as xr

    N = 1024 # number of grid points
    step = 300 # time step size [s]

    xmin = 0.
    xmax = 28e6
    xdel = xmax - xmin
    # Grid space coordinates
    x = np.linspace(xmin, xmax, N, endpoint=False)
    lon = x / xdel * 360 - 180
    # Spectral space wavenumbers
    k = 2 * np.pi * np.arange(0, N//2+1, dtype=complex) / xdel

    tau    = 10. * DAYS # relaxation time scale [s]
    kappa  = 3.26e5 # diffusion coefficient [m²/s]
    linear = - 1. / tau - kappa * k**2

    # Stationary LWA
    A0 = 19 * hann(lon, 10., 135)
    # Uniform constant forcing throughout the domain
    base_forcing = 1.825e-5 * np.ones_like(x)

    def flux(lwa):
        return (args.urefcg - 2. * args.alpha * A0 - args.alpha * lwa) * lwa

    # Model with base forcing only for spin-up
    def nonlinear(lwa_spectral, t):
        lwa = np.fft.irfft(lwa_spectral)
        fluxdiv = 1j * k * np.fft.rfft(flux(lwa))
        forcing = base_forcing
        return np.fft.rfft(forcing) - fluxdiv
    model = ETDRK4Diagonal(step, linear, nonlinear)
    # Start with no transient LWA
    init = np.fft.rfft(np.zeros_like(A0))
    # Run 200 days of spin-up to reach steady state
    _, [_, init] = model.run(init, 200*DAYS//step, t=0, keep_every=0)

    lwa, flx = [], []
    # Run each ensemble members with a different forcing strength
    strengths = np.linspace(0.0000, 0.0007, args.nens)
    for strength in strengths:
        # Model with upstream forcing
        def nonlinear(lwa_spectral, t):
            lwa = np.fft.irfft(lwa_spectral)
            fluxdiv = 1j * k * np.fft.rfft(flux(lwa))
            forcing = base_forcing + (strength * hann(t, args.forcing_peak*DAYS, 3.5*DAYS) * hann(lon, -75., 75.))
            return np.fft.rfft(forcing) - fluxdiv
        model = ETDRK4Diagonal(step, linear, nonlinear)
        # Run 8 days of actual simulation
        t, A = model.run(init, 8*DAYS//step, t=-5*DAYS, keep_every=72)
        # Add simulation data to ensemble
        A = np.fft.irfft(A)
        lwa.append(A + A0)
        flx.append(flux(A))

    data_vars = {
        "lwa": (("number", "time", "longitude"), np.asarray(lwa, dtype=np.float32)),
        "flux": (("number", "time", "longitude"), np.asarray(flx, dtype=np.float32)),
        "lwa_stat": (("longitude",), A0),
        "lwa_thresh": (("longitude",), 0.5 * (args.urefcg - 2. * args.alpha * A0) / args.alpha)
    }
    coords = {
        "number": np.arange(len(lwa)),
        "time": t / DAYS,
        "longitude": lon,
    }
    attrs = {
        "nens": args.nens,
        "urefcg": args.urefcg,
        "alpha": args.alpha,
        "forcing_peak": args.forcing_peak,
        "step": step,
        "tau": tau,
        "kappa": kappa,
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    }
    xr.Dataset(data_vars, coords, attrs).to_netcdf(args.outfile)

