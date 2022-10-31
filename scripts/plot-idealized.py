import argparse 
import json

parser = argparse.ArgumentParser()
parser.add_argument("infiles", type=str, nargs="+")
parser.add_argument("--output", type=str, default=None)
parser.add_argument("--data-output", type=str, default=None)

if __name__ == "__main__":
    args = parser.parse_args()

    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    from .common import plotting, esa

    # Upstream region
    up_d, up_s, up_e = -2, -90, -15
    # Blocking region
    bl_d, bl_s, bl_e =  0, -30,  10

    def prep_data(ds):
        # 2° resolution is good enough for these plots
        ds = ds.interp({ "longitude": np.arange(-104, 77, 2) })
        # The prescribed stationary wave activity leads to a stationary component
        # of the transient LWA flux. Remove this stationary component by
        # subtracting the stationary state after spin-up at the beginning of the
        # simulation from the flux field, then visualize flux anomalies.
        ds["flux*"] = ds["flux"] - ds["flux"].isel(time=0)
        # Determine if member is jammed at block time
        lwa_trans = ds["lwa"] - ds["lwa_stat"]
        ds["jammed0"] = (lwa_trans.sel({ "time": bl_d              }) > ds["lwa_thresh"]).any(dim="longitude")
        ds["jammed1"] = (lwa_trans.sel({ "time": slice(bl_d, None) }) > ds["lwa_thresh"]).any(dim=["longitude", "time"])
        return ds

    # Mean of the reference ensmeble
    mean = xr.load_dataset(args.infiles[0])
    mean = prep_data(mean)
    nens = mean["number"].size
    mean = mean.mean(dim="number")
    # Full ensemble including perturbed members
    data = []
    labels = []
    for infile in args.infiles:
        _ds = xr.load_dataset(infile)
        if not labels:
            labels.append("Reference")
            ref_alpha = _ds.attrs["alpha"]
            ref_urefcg = _ds.attrs["urefcg"]
            ref_peak = _ds.attrs["forcing_peak"]
        elif _ds.attrs["alpha"] != ref_alpha:
            labels.append(r"$\alpha = {:.2f}$".format(_ds.attrs["alpha"]))
        elif _ds.attrs["urefcg"] != ref_urefcg:
            labels.append(r"$u_\mathrm{{ref}} + c_\mathrm{{g}} = {:.0f} \;\mathrm{{m}}\,\mathrm{{s}}^{{-1}}$".format(_ds.attrs["urefcg"]))
        elif _ds.attrs["forcing_peak"] != ref_peak:
            labels.append(r"$t_\mathrm{{up}} = {:.1f} \;\mathrm{{d}}$".format(_ds.attrs["forcing_peak"]))
        else:
            labels.append("???")
        data.append(prep_data(_ds))
    data = xr.concat(data, dim="number")
    # Aggregated data for the sensitivity analysis
    upstream = data.sel({ "time": up_d, "longitude": slice(up_s, up_e) }).mean(dim="longitude")
    blocked  = data.sel({ "time": bl_d, "longitude": slice(bl_s, bl_e) }).mean(dim="longitude")

    fig = plt.figure(figsize=(9, 4.1))

    lon  = data["longitude"].values
    time = data["time"].values

    _, ax_lwa = plotting.add_hov_axes(fig, (0.02, 0.225, 0.16, 0.71), lonlim=(-105, 75), latlim=(35, 65))
    _, ax_flx = plotting.add_hov_axes(fig, (0.20, 0.225, 0.16, 0.71), lonlim=(-105, 75), latlim=(35, 65))
    _, ax_esa = plotting.add_hov_axes(fig, (0.38, 0.225, 0.16, 0.71), lonlim=(-105, 75), latlim=(35, 65))
    ax_sct = fig.add_axes((0.67, 0.115, 0.31, 0.82))
    cx_lwa = fig.add_axes((0.02, 0.12, 0.16, 0.025))
    cx_flx = fig.add_axes((0.20, 0.12, 0.16, 0.025))
    cx_esa = fig.add_axes((0.38, 0.12, 0.16, 0.025))

    ct_lwa = ax_lwa.contour(lon, time, mean["lwa"].values, levels=[50], colors="#000", linewidths=0.8)
    plt.clabel(ct_lwa, fontsize="small")
    cf_lwa = ax_lwa.contourf(
        lon, time, mean["lwa"].values,
        levels=[0, 10, 20, 30, 40, 50, 60, 70],
        colors=plotting.lwa_colors,
        extend="max"
    )
    ax_lwa.plot([bl_s, bl_e], [bl_d, bl_d], color="#000", linewidth=2)
    cb_lwa = fig.colorbar(cf_lwa, cax=cx_lwa, orientation="horizontal", spacing="proportional", extendfrac=0.1)
    cb_lwa.set_ticks([0, 20, 40, 60])
    cb_lwa.set_label(r"$A_0 + \hat A$ [$\mathrm{m} \; \mathrm{s}^{-1}$]")
    ax_lwa.set_yticklabels([])

    ct_flx = ax_flx.contour(lon, time, mean["flux*"].values, levels=[360], colors="#000", linewidths=0.8)
    plt.clabel(ct_flx, fontsize="small")
    cf_flx = ax_flx.contourf(
        lon, time, mean["flux*"].values,
        levels=[120, 240, 360, 480],
        colors=plotting.flux_colors[4:],
        extend="both"
    )
    ax_flx.plot([up_s, up_e], [up_d, up_d], color="#000", linewidth=2)
    cb_flx = fig.colorbar(cf_flx, cax=cx_flx, orientation="horizontal", spacing="proportional", extendfrac=0.1)
    cb_flx.set_ticks([120, 240, 360, 480])
    cb_flx.set_label(r"$\hat F$ anomaly [$\mathrm{m}^2 \; \mathrm{s}^{-2}$]")
    ax_flx.set_yticklabels([])

    ct_esa = ax_esa.contour(lon, time, mean["flux*"].values, levels=[360], colors="#000", linewidths=0.8)
    plt.clabel(ct_esa, fontsize="small")
    target = blocked["lwa"].values
    source = data["flux"].values
    cf_esa = ax_esa.contourf(
        lon, time, esa.esa_corr(source, target),
        levels=plotting.esa_levels,
        colors=plotting.esa_colors,
        extend="both"
    )
    cb_esa = fig.colorbar(cf_esa, cax=cx_esa, orientation="horizontal", spacing="proportional", extendfrac=0.1)
    cb_esa.set_ticks([-0.8, -0.4, 0.4, 0.8])
    cb_esa.set_label(r"Sensitivity to $\hatF$")
    ax_esa.yaxis.set_major_formatter(r"${x:+.0f} \;\mathrm{{d}}$")

    _colors = ["#999", "#999", "#CCC", "#CCC"]
    for i in range(data["number"].size // nens):
        slc = slice(i*nens, (i+1)*nens)
        src = upstream["flux*"].values[slc]
        jam0 = data["jammed0"].values[slc]
        jam1 = data["jammed1"].values[slc]
        tgt = target[slc]
        # Reference simulation
        if i == 0:
            kw = { "color": "#000", "zorder": 10 }
            fc = np.full(jam1.shape, "#FFF")
            fc[jam1] = "#CCC"
            fc[jam0] = "#000"
            ax_sct.scatter(src, tgt, s=30, marker="o", color="#000", facecolor=fc, zorder=11)
            # Save values of src-tgt curve to file for use in other plots
            curve = { "x": src.tolist(), "y": tgt.tolist(), "z0": jam0.tolist(), "z1": jam1.tolist() }
            if args.data_output is None:
                print(json.dumps(curve, indent=4))
            else:
                with open(args.data_output, "w") as f:
                    json.dump(curve, f, indent=4)
        # Perturbation experiments
        else:
            kw = {
                "color": _colors.pop(),
                "zorder": 0,
                "label": labels[i],
                "linestyle": ("solid" if i % 2 == 0 else "dashed")
            }
        kw["zorder"] -= 1
        ax_sct.plot(src, tgt, linewidth=1.2, **kw)

    handles = [
        mlines.Line2D([], [], marker="o", color="#000", linewidth=0, label=r"$\hat A > {\hat A}_\mathrm{C}$ ($+0$ d)", markerfacecolor="#000"),
        mlines.Line2D([], [], marker="o", color="#000", linewidth=0, label=r"$\hat A > {\hat A}_\mathrm{C}$", markerfacecolor="#CCC"),
        mlines.Line2D([], [], marker="o", color="#000", linewidth=0, label=r"$\hat A < {\hat A}_\mathrm{C}$", markerfacecolor="#FFF"),
    ]
    leg1 = ax_sct.legend(handles=handles, loc="upper left", handlelength=1.2, fontsize="small")
    leg2 = ax_sct.legend(loc="lower center", fontsize="small", handlelength=1.2, ncol=2)
    ax_sct.add_artist(leg1)

    ax_sct.set_xlim((-50, 950))
    ax_sct.set_xlabel(r"agg. LWA Flux (${:.0f} \;\mathrm{{d}}$) [$\mathrm{{m}}^2 \; \mathrm{{s}}^{{-2}}$]".format(up_d))
    ax_sct.set_ylim((22, 84))
    ax_sct.set_ylabel(r"agg. LWA (${:.0f} \;\mathrm{{d}}$) [$\mathrm{{m}} \; \mathrm{{s}}^{{-1}}$]".format(bl_d))
    plotting.set_grid(ax_sct)

    ax_lwa.set_title("(a) LWA", loc="left")
    ax_flx.set_title("(b) LWA Flux", loc="left")
    ax_esa.set_title("(c) ESA (corr.)", loc="left")
    ax_sct.set_title("(d) Sensitivity Scatter", loc="left")

    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output, dpi=200)


