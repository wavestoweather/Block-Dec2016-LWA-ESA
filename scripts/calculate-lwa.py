import argparse
import os
import re
import datetime as dt

import numpy as np

from .common.regrid import half_resolution


def pretty_dt(dt):
    seconds = dt.total_seconds()
    minutes = int(seconds / 60)
    seconds = int(seconds - 60 * minutes)
    return "{}:{:0>2}".format(minutes, seconds)


def encoding_int16(scale_factor, add_offset):
    return {
        "dtype": "int16",
        "scale_factor": scale_factor,
        "add_offset": add_offset,
        "_FillValue": -32768,
        "zlib": True
    }


class QuasiGeostrophicCalculator:

    # Attributes for the netcdf output files
    attrs = {
        "framework": "QG"
    }
    # File name suffix
    suffix = "-qg"

    def __init__(self, budget=False):
        self.budget = budget
        # Variables (in order) returned by calculate()
        self.return_vars = ["pv", "lwa", "f1", "f2", "f3"]
        # Encoding information for netcdf output
        self.encoding = {
            "lwa": encoding_int16(0.02, 400),
            "pv": encoding_int16(7.0e-8, 0),
            "f1": encoding_int16(0.15, 0),
            "f2": encoding_int16(0.15, 0),
            "f3": encoding_int16(0.15, 0),
        }
        # More fields returned if budget is included
        if budget:
            self.return_vars.extend(["czaf", "demf", "mhf"])
            self.encoding["czaf"] = encoding_int16(1.0e-6, 0.)
            self.encoding["demf"] = encoding_int16(1.0e-6, 0.)
            self.encoding["mhf"] =  encoding_int16(1.0e-6, 0.)

    def calculate(self, xlon, ylat, plev, u, v, t):
        # On-demand import
        from .hn2016_falwa.oopinterface import QGField
        nlat = ylat.size
        # QGField has specific requirements for ordering of the dimensions,
        # need to flip some to make the data fit
        assert ylat[0] > ylat[1]
        assert plev[0] < plev[1]
        # Need to flip sign of latitude and meridional wind on mirrored
        # hemisphere to obtain same LWA, zonal wind and temperature stay
        # unchanged. The impact of the discontinuity created due to the flipped
        # sign in the meridional wind at the equator does not have much impact
        # on LWA or the mid-latitude fluxes.
        plev = plev[::-1]
        ylat = np.concatenate([-ylat, ylat[-2::-1]])
        v    = np.concatenate([-v[::-1,:,:], v[::-1,-2::-1,:]], axis=1)
        u    = np.concatenate([ u[::-1,:,:], u[::-1,-2::-1,:]], axis=1)
        t    = np.concatenate([ t[::-1,:,:], t[::-1,-2::-1,:]], axis=1)
        # Up to 25 km = 28 hPa < 10 hPa, 1 km level spacing, 7 km scale height
        kmax   = 26
        dz     = 1000.
        hscale = 7000.
        # Calculate prefactor (see commit 2de8e284 and issue #39 in hn2016_falwa)
        prefactor = np.sum(dz * np.exp(-np.arange(1, kmax-1) * dz / hscale))
        # Everything ready for QGField to do the work
        fld = QGField(xlon, ylat, plev, u, v, t, scale_height=hscale, dz=dz, kmax=kmax, prefactor=prefactor)
        fld.interpolate_fields()
        fld.compute_reference_states(northern_hemisphere_results_only=True)
        fld.compute_lwa_and_barotropic_fluxes(northern_hemisphere_results_only=True)
        # Flip fields back to original orientation
        out = (
            -fld.qgpv[7,:nlat,:],
            fld.lwa_baro[::-1,:],
            fld.adv_flux_f1[::-1,:],
            fld.adv_flux_f2[::-1,:],
            fld.adv_flux_f3[::-1,:]
        )
        # Add budget terms to output if requested
        if self.budget:
            out = out + (
                fld.convergence_zonal_advective_flux[::-1,:],
                fld.divergence_eddy_momentum_flux[::-1,:],
                fld.meridional_heat_flux[::-1,:],
            )
        return out

    def get_var_attrs(self, var):
        if var == "pv":
            return { "name": "QGPV", "unit": "1/s", "level": "7 km" }
        elif var == "lwa":
            return { "name": "<A>cos(φ)", "unit": "m/s", }
        elif var in ["f1", "f2", "f3"]:
            return { "name": var.upper(), "unit": "m²/s²" }
        elif var == "czaf":
            return { "name": "Convergence of the Zonal Advective Flux", "unit": "m/s²" }
        elif var == "demf":
            return { "name": "Divergence of the Eddy Momentum Flux", "unit": "m/s²" }
        elif var == "mhf":
            return { "name": "Meridional Heat Flux", "unit": "m/s²" }
        else:
            return {}



INIT_REGEX = re.compile(r".*EVAL-([1-9][0-9][0-9][0-9]-[0-9][0-2]-[0-3][0-9]T[0-2][0-9])Z.nc")
# In case of forecast evaluation, time in the dataset (normally cotaining valid
# time) is replaced with forecast initialization time (this way processing
# below falls back to the scheme for reanalysis data and no additional
# exception is needed). The valid time is kept in a new coordinate.
def preprocess_eval(ds):
    filename = ds.encoding["source"] # https://github.com/pydata/xarray/issues/2550
    match = INIT_REGEX.match(filename)
    assert match is not None
    init = dt.datetime.strptime(match.groups()[0], "%Y-%m-%dT%H")
    return ds.assign_coords({ "time": [init], "valid": ds["time"].values })


parser = argparse.ArgumentParser()
parser.add_argument("--eval", action="store_true", help="""
    Process in forecast evaluation mode (merges multiple forecasts for the same
    valid time)
""")
parser.add_argument("--budget", action="store_true", help="""
    Store budget terms in output (only available in QG framework).
""")
parser.add_argument("--half-resolution", action="store_true", help="""
    Regrid the input data in the horizonatl to a mesh with half the resolution.
""")
parser.add_argument("infile", nargs="+", type=str, help="""
    Input files (u, v, t) for the computation.
""")

if __name__ == "__main__":
    
    args = parser.parse_args()

    # Import big modules after argument parsing so --help is fast
    import pandas as pd
    import xarray as xr

    calculator = QuasiGeostrophicCalculator(args.budget)

    mf_kwargs = {
        "parallel": False,
        "chunks": { "time": 1 },
    }
    # Combine multiple forecasts for the same valid time for evaluation
    if args.eval:
        mf_kwargs["combine"] = "nested"
        mf_kwargs["concat_dim"] = "time"
        mf_kwargs["preprocess"] = preprocess_eval
    # Open dataset and perform basic integrity check
    data = xr.open_mfdataset(args.infile, **mf_kwargs)
    assert "u" in data, "zonal wind data not found in input files"
    assert "v" in data, "meridional wind data not found in input files"
    assert "t" in data, "temperature data not found in input files"

    # Boolean indicating if this is ensemble data
    is_ens = "number" in data

    # Pressure levels in hPa
    p = data.level.values
    # Latitude and longitude grid
    lat = data.latitude.values
    lon = data.longitude.values

    init = pd.to_datetime(data.time.values)[0]
    start = init + dt.timedelta(hours=48)
    # Remove first 48 hours of forecasts (not reliable for ESA) but keep for
    # forecast evaluation purposes
    if is_ens and not args.eval:
        data = data.sel(time=slice(start, None))

    out = { var: [] for var in calculator.return_vars }

    # Take time for progress reports during the computation
    _t_start = dt.datetime.now()

    valids = pd.to_datetime(data.time.values)
    assert len(valids) > 0, valids
    for i, time in enumerate(valids):

        # Report on progress
        _elapsed = dt.datetime.now() - _t_start
        print("{init:%Y-%m-%dT%H}Z +{lead:0>3} (step {step}/{total}, {elapsed} elapsed, {remaining} remaining)".format(
            init=init,
            lead=int((time-init).total_seconds() / 3600),
            step=i+1,
            total=len(valids),
            elapsed=pretty_dt(_elapsed),
            remaining=pretty_dt((_elapsed / i) * (len(valids) - i)) if i != 0 else "?:??"
        ))

        if is_ens:
            for number in data.number.values:
                u = data.u.sel(time=time, number=number).values
                v = data.v.sel(time=time, number=number).values
                t = data.t.sel(time=time, number=number).values
                # Reduce the resolution of the input data if desired. Make
                # sure to call with coordinates of original grid, lat and lon
                # is overwritten after the first iteration
                if args.half_resolution:
                    lat, lon, u, v, t = half_resolution(data.latitude.values, data.longitude.values, u, v, t)
                # Compute output fields
                fields = calculator.calculate(lon, lat, p, u, v, t)
                # Add each variable in the output fields to its collection
                for var, field in zip(calculator.return_vars, fields):
                    assert field.shape == (lat.size, lon.size)
                    out[var].append(field)
        else:
            u = data.u.sel(time=time).values
            v = data.v.sel(time=time).values
            t = data.t.sel(time=time).values
            if args.half_resolution:
                lat, lon, u, v, t = half_resolution(data.latitude.values, data.longitude.values, u, v, t)
            fields = calculator.calculate(lon, lat, p, u, v, t)
            for var, field in zip(calculator.return_vars, fields):
                assert field.shape == (lat.size, lon.size)
                out[var].append(field)

    print("{init:%Y-%m-%dT%H}Z writing...".format(init=init))

    # Reshape flattened output fields to appropriate size
    outshape = (data.time.size,) + ((data.number.size,) if is_ens else tuple()) + (lat.size, lon.size)
    # The corresponding dimension names
    dims = ("time",) + (("number",) if is_ens else tuple()) + ("latitude", "longitude")
    # Assemble output arrays for writing
    data_vars = {
        var: (dims, np.asarray(out[var]).reshape(outshape))
        for var in calculator.return_vars
    }

    # Coordinates of output fields depend on resolution choice
    coords = {
        "time": data.time,
        "latitude": data.latitude[::2] if args.half_resolution else data.latitude,
        "longitude": data.longitude[::2] if args.half_resolution else data.longitude
    }
    # Make sure coordinates match what those from the processing (lat, lon
    # variables not used directly, since data.latitude and data.longitude carry
    # metadata, which should be preserved)
    assert np.allclose(coords["latitude"].values, lat)
    assert np.allclose(coords["longitude"].values, lon)
    # Need to add ensemble coordinates if input data is from ENS
    if is_ens:
        coords["number"] = data.number
    if args.eval:
        coords["valid"] = data.valid

    # Global attributes
    attrs = {
        "created": dt.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    }
    attrs.update(calculator.attrs)

    # Combine fields, coordinates and general metadata into Dataset, then set
    # variable-specific metadata
    out = xr.Dataset(data_vars, coords, calculator.attrs)
    for var in calculator.return_vars:
        for key, value in calculator.get_var_attrs(var).items():
            out[var].attrs[key] = value

    # Generate output file name
    name = "data/"
    if args.eval:
        name += "EVAL-{:%Y-%m-%dT%H}Z".format(pd.to_datetime(data.valid.values[0]))
    elif is_ens:
        name += "ENS-{:%Y-%m-%dT%H}Z".format(init)
    else:
        name += "ERA-{:%Y-%m-%d}-to-{:%Y-%m-%d}".format(
            pd.to_datetime(data.time.values[0]),
            pd.to_datetime(data.time.values[-1]),
        )
    name += "-lwa" + calculator.suffix + ".nc"

    # Write to disk
    out.to_netcdf(name, encoding=calculator.encoding)

    print("{init:%Y-%m-%dT%H}Z finished in {elapsed}!".format(
        init=init,
        elapsed=pretty_dt(dt.datetime.now() - _t_start)
    ))

