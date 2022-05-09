import argparse
import json
from .plot import verify_path, parse_target, parse_datehour


parser = argparse.ArgumentParser()
parser.add_argument("reanalysis", type=verify_path, help="""
    Path to reanalysis dataset.
""")
parser.add_argument("start", type=parse_datehour, help="""
    Start date of integration period.
""")
parser.add_argument("target", type=parse_target, help="""
    End date of integration period and target box for aggregation.
""")
parser.add_argument("--output", type=str, default=None, help="""
    Output file path. If omitted plt.show is called instead.
""")
parser.add_argument("--statistics-output", type=str, default=None, help="""
    Output file path for statistical data collected during plotting. If
    omitted, data is printed to stdout.
""")

if __name__ == "__main__":
    args = parser.parse_args()

    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    from .common.plotting import add_map_axes, plot_box, transform, format_date_full
    from .common.texify import texify

    ana = xr.open_dataset(args.reanalysis)
    assert "lwa" in ana, "no LWA fields in dataset"
    assert "czaf" in ana, "no convergence of the zonal advective flux fields in dataset"
    assert "demf" in ana, "no divergence of the eddy momentum flux fields in dataset"
    assert "mhf" in ana, "no meridional heat flux fields in dataset"

    # Only keep data from evaluation period
    ana = ana.sel({ "time": slice(args.start, args.target.valid) })

    # Determine LWA delta in given period
    lwa_start = ana["lwa"].sel({ "time": args.start })
    lwa_end   = ana["lwa"].sel({ "time": args.target.valid })
    dlwa = lwa_end - lwa_start
    # Integrate budget terms over given time period
    quad = ana.integrate(coord="time", datetime_unit="s")
    # Obtain residual from unexplained LWA change in budget
    quad["lwa"] = dlwa
    quad["res"] = quad["lwa"] - quad["czaf"] - quad["demf"] - quad["mhf"]

    target_slice = {
        "latitude": slice(args.target.N, args.target.S),
        "longitude": slice(args.target.W, args.target.E)
    }
    def stat_target(dterm, dlwa):
        dt = dterm.sel(target_slice)
        dl = dlwa.sel(target_slice)
        corr = np.corrcoef(dt.values.flatten(), dl.values.flatten())[0,1]
        agg  = dt.weighted(np.cos(np.deg2rad(dt["latitude"]))).mean(dim=["latitude", "longitude"])
        return { "corr": corr, "delta": float(agg) }

    stats = { var: stat_target(quad[var], quad["lwa"]) for var in quad }
    if args.statistics_output is None:
        print(json.dumps(stats, indent=4))
    else:
        with open(args.statistics_output, "w") as f:
            json.dump(stats, f, indent=4)

    fig = plt.figure(figsize=(4, 5.12))
    ax_lwa  = add_map_axes(fig, (0.14, 0.75, 0.84, 0.23), noxlabels=True)
    ax_czaf = add_map_axes(fig, (0.14, 0.47, 0.84, 0.23), noxlabels=True)
    ax_res  = add_map_axes(fig, (0.14, 0.19, 0.84, 0.23))
    cx_lwa  = fig.add_axes((0.14, 0.095, 0.84, 0.023))

    lon = quad["longitude"].values
    lat = quad["latitude"].values

    dlwa_kwargs = {
        "colors": ['#01665e', '#35978f', '#80cdc1', '#c7eae5', '#f5f5f5', '#f6e8c3', '#dfc27d', '#bf812d', '#8c510a'],
        "levels": [-120, -80, -40, -20, 20, 40, 80, 120],
        "extend": "both",
        "transform": transform,
    }
    flux_kwargs = {
        "levels": [-2000, -1600, -1200, -800, -400, 400, 800, 1200, 1600, 2000],
        "linewidths": 0.8,
        "colors": "#000",
    }

    ax_lwa.set_title("(a) $\Delta$" + texify(ana["lwa"].attrs["name"]), loc="left")
    cf = ax_lwa.contourf(lon, lat, quad["lwa"], **dlwa_kwargs)

    cb = fig.colorbar(cf, cax=cx_lwa, spacing="proportional", orientation="horizontal")
    cb.set_label(r"LWA Change [" + texify("m/s") + "]")
    cb.set_ticks([-120, -80, -40, 0, 40, 80, 120])

    _start = format_date_full(args.start).replace("2016 ", "")
    _end   = format_date_full(args.target.valid).replace("2016 ", "")
    ax_lwa.set_title("{} to {}".format(_start, _end), loc="right", fontsize="small")

    _flambda = texify("F1+F2+F3")
    ax_czaf.set_title("(b) " + _flambda + ", " + _flambda + " Convergence", loc="left")
    ax_czaf.set_title("$r = {corr:.2f}$".format(**stats["czaf"]), loc="right", fontsize="small")
    ax_czaf.contourf(lon, lat, quad["czaf"], **dlwa_kwargs)
    ct = ax_czaf.contour(lon, lat, (ana["f1"] + ana["f2"] + ana["f3"]).mean(dim="time"), **flux_kwargs)
    plt.clabel(ct, fontsize="x-small")

    ax_res.set_title("(c) Non-Conservative", loc="left")
    ax_res.set_title("$r = {corr:.2f}$".format(**stats["res"]), loc="right", fontsize="small")
    ax_res.contourf(lon, lat, quad["res"], **dlwa_kwargs)

    for ax in [ax_lwa, ax_czaf, ax_res]:
        plot_box(ax, args.target, linewidth=2)

    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output, dpi=200)

