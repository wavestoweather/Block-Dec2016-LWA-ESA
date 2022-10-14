import argparse
import collections
import datetime as dt
import operator
import os
import json
from functools import partial, reduce


def parse_date(s):
    return dt.datetime.strptime(s, "%Y-%m-%d")

def parse_datehour(s):
    return dt.datetime.strptime(s, "%Y-%m-%dT%H")

def verify_path(p):
    assert os.path.exists(p), "'{}' does not exist".format(p)
    return p

TargetProperties = collections.namedtuple("TargetProperties", ["valid", "N", "E", "S", "W"])
def parse_target(s):
    valid, *trbl = s.split(",")
    t, r, b, l = [float(x) for x in trbl]
    if t < b:
        t, b = b, t
    if r < l:
        r, l = l, r
    return TargetProperties(parse_datehour(valid), t, r, b, l)

def parse_metric(s):
    assert s == "average" or s.startswith("area>")
    return s

def parse_delta_hours(s):
    return [int(dh) for dh in s.split(",")]

parser = argparse.ArgumentParser()
# Optional arguments
parser.add_argument("--reanalysis", type=verify_path, default=None, help="""
    Point this at a reanalysis dataset to produce an additional panel with
    reanalysis LWA and PV as well as a reanalysis timeseries in the plume panel
""")
parser.add_argument("--ncluster", type=int, default=5, help="""
    How many members selected for the cluster analysis
""")
parser.add_argument("--nsamples", type=int, default=1000, help="""
    Number of samples drawn for the bootstrap significance test
""")
parser.add_argument("--sample-size", type=int, default=None, help="""
    ...
""")
parser.add_argument("--significance-level", type=float, default=0.1, help="""
    Threshold for standard deviation of correlation distribution below which
    results are declared as statistically significant.
""")
parser.add_argument("--target-metric", type=parse_metric, default="average", help="""
    The target metric definition. Available choices: 'average' (area-weighted
    average of LWA inside the target region, default), `area>X' (size of region
    where LWA exceeds X inside the target region).
""")
parser.add_argument("--source", type=str, default="lwa", help="""
    The source variable for the sensitivity analysis. Affects the sensitivity
    plots in map(s), hovmoeller and scatter.
""")
parser.add_argument("--scatter", type=parse_target, default=None, help="""
    valid,N,E,S,W for scatter plot source definition.
""")
parser.add_argument("--scatter-idealized", type=verify_path, default=None, help="""
    Reference curve added to scatter plot (json format).
""")
parser.add_argument("--hovmoeller-extent", type=str, default="target", help="""
    The meridional extend used in the spatial reduction of data for the
    Hovmöller diagram. Specify as N,S or use the special keywords 'target'
    (default) or 'scatter' to use the extend of the target or scatter regions,
    respectively.
""")
parser.add_argument("--delta-hours", type=parse_delta_hours, default="0,-24,-48,-72,-96", help="""
    Which valid time deltas to show in the maps
""")
parser.add_argument("--panel-letters", default=None, help="""
    Override and force panel labels
""")
parser.add_argument("--start", type=parse_datehour, default=None, help="""
    Restrict the plotted period at the start
""")
parser.add_argument("--end", type=parse_datehour, default=None, help="""
    Restrict the plotted period at the end
""")
parser.add_argument("--resample", type=str, default=None, help="""
    Resample the dataset. Accepts resample indexers for xarray.Dataset.resample
""")
parser.add_argument("--output", type=str, default=None, help="""
    Output file path. If omitted plt.show is called instead.
""")
# Mandatory arguments
parser.add_argument("layout", type=str, help="""
    Choose the contents and layout of the figure
""")
parser.add_argument("target", type=parse_target, help="""
    Target definition: valid,N,E,S,W. Valid date in YYYY-MM-DDTHH format, box
    definition in degrees latitude/longitude, values in CSS order).
""")
parser.add_argument("data", type=verify_path, nargs="+", help="""
    Ensemble data sets that are combined
""")


if __name__ == "__main__":

    args = parser.parse_args()

    # Delayed import of big modules so --help is fast
    import numpy as np
    import xarray as xr
    import matplotlib.pyplot as plt
    from scipy.stats import linregress

    from .common import metrics, plotting
    from .common.texify import texify
    from .common.plotting import (
        add_map_axes, add_hov_axes, plot_box, set_grid, date_formatter,
        format_date_full, mult_limit, esa_blue, esa_red
    )

    # Panel naming: a), b), ...
    letters = list(reversed(
        "abcdefghijklmnop" if args.panel_letters is None else args.panel_letters
    ))
    def set_panel_title(ax, title, **kwargs):
        try:
            ax.set_title("({}) {}".format(letters.pop(), title), loc="left", **kwargs)
        except IndexError:
            ax.set_title(title, loc="left", **kwargs)


    # Time deltas for map panels can be specified by the user, but may be
    # overwritten by specific plot presets
    delta_hours = args.delta_hours

    # Meridional extent for Hovmöller data aggregation (required here already
    # to for the axes layout in add_hov_axes). Data is given with latitude
    # coordinates from N to S, so specify slice in same order.
    if args.hovmoeller_extent == "target":
        hovmoeller_extent = slice(args.target.N, args.target.S)
    # or take from definition of scatter aggregation region...
    elif args.hovmoeller_extent == "scatter":
        assert args.scatter is not None, "scatter region is not specified for use in hovmoeller reduction"
        hovmoeller_extent = slice(args.scatter.N, args.scatter.S)
    # or parse directly given extent
    else:
        _hovNS = [float(x) for x in args.hovmoeller_extent.split(",")]
        assert len(_hovNS) == 2, "Meridional extent for Hovmöller data reduction must be specified as N,S"
        hovmoeller_extent = slice(max(_hovNS), min(_hovNS))

    # The panel layout is created first. Panels are given specific names, which
    # are re-used in different layouts, so they can be filled by the same
    # routines below.

    # 3-panel plot with target plume and LWA/PV rank clusters
    if args.layout == "plume":
        fig = plt.figure(figsize=(9, 3.5))
        ax_plu = fig.add_axes(     (0.06, 0.08, 0.45, 0.83))
        ax_top = add_map_axes(fig, (0.58, 0.58, 0.35, 0.33), lonlim=(-105, 75))
        ax_bot = add_map_axes(fig, (0.58, 0.08, 0.35, 0.33), lonlim=(-105, 75))
        cx_lwa = fig.add_axes((0.925, 0.08, 0.01, 0.83)) # for the LWA colorbar
    # Plume/cluster panels with reanalysis LWA/PV map on top
    elif args.layout == "plume+reanalysis":
        fig = plt.figure(figsize=(9, 6.5))
        ax_ana = add_map_axes(fig, (0.06, 0.57, 0.83, 0.38), lonlim=(-135, 75))
        ax_plu = fig.add_axes(     (0.06, 0.05, 0.51, 0.42))
        ax_top = add_map_axes(fig, (0.65, 0.29, 0.34, 0.18), lonlim=(-105, 75), noxlabels=True)
        ax_bot = add_map_axes(fig, (0.65, 0.05, 0.34, 0.18), lonlim=(-105, 75))
        cx_lwa = fig.add_axes((0.91, 0.57, 0.013, 0.38)) # for the LWA colorbar
    # Plume/cluster panels with Hovmöller ESA plot
    elif args.layout == "plume+hovmoeller":
        fig = plt.figure(figsize=(9, 6.5))
        ax_plu = fig.add_axes(     (0.06, 0.65, 0.42, 0.30))
        ax_top = add_map_axes(fig, (0.06, 0.34, 0.38, 0.20), lonlim=(-105, 75))
        ax_bot = add_map_axes(fig, (0.06, 0.05, 0.38, 0.20), lonlim=(-105, 75))
        _, ax_hov = add_hov_axes(
            fig, (0.54, 0.16, 0.38, 0.79),
            latlim=(hovmoeller_extent.stop, hovmoeller_extent.start), # call with S first
            relref=args.target.valid
        )
        # Colorbars for LWA and correlation
        cx_lwa = fig.add_axes((0.45, 0.05, 0.013, 0.49))
        cx_esa = fig.add_axes((0.54, 0.08, 0.38, 0.020))
    # Plume/cluster panels with 4 ESA maps
    elif args.layout == "plume+maps":
        fig = plt.figure(figsize=(9, 6.5))
        ax_plu = fig.add_axes(     (0.06, 0.65, 0.42, 0.30))
        ax_top = add_map_axes(fig, (0.06, 0.34, 0.38, 0.20), lonlim=(-105, 75))
        ax_bot = add_map_axes(fig, (0.06, 0.05, 0.38, 0.20), lonlim=(-105, 75))
        ax_maps = [
            add_map_axes(fig, (0.61, 0.80, 0.37, 0.15), noxlabels=True),
            add_map_axes(fig, (0.61, 0.59, 0.37, 0.15), noxlabels=True),
            add_map_axes(fig, (0.61, 0.38, 0.37, 0.15), noxlabels=True),
            add_map_axes(fig, (0.61, 0.17, 0.37, 0.15))
        ]
        # Colorbars for LWA and correlation
        cx_lwa = fig.add_axes((0.45, 0.05, 0.013, 0.49))
        cx_esa = fig.add_axes((0.61, 0.08, 0.37, 0.020))
    # 4 ESA maps and Hovmöller
    elif args.layout == "maps+hovmoeller":
        fig = plt.figure(figsize=(9, 6.5))
        ax_maps = [
            add_map_axes(fig, (0.06, 0.78, 0.45, 0.17), noxlabels=True),
            add_map_axes(fig, (0.06, 0.54, 0.45, 0.17), noxlabels=True),
            add_map_axes(fig, (0.06, 0.30, 0.45, 0.17), noxlabels=True),
            add_map_axes(fig, (0.06, 0.06, 0.45, 0.17))
        ]
        _, ax_hov = add_hov_axes(
            fig, (0.54, 0.16, 0.38, 0.79),
            latlim=(hovmoeller_extent.stop, hovmoeller_extent.start), # call with S first
            relref=args.target.valid
        )
        cx_esa = fig.add_axes((0.54, 0.08, 0.38, 0.020))
    # Column of 5 ESA maps
    elif args.layout == "maps" or args.layout == "maps5":
        fig = plt.figure(figsize=(4, 8))
        ax_maps = [
            add_map_axes(fig, (0.14, 0.84, 0.84, 0.15), noxlabels=True),
            add_map_axes(fig, (0.14, 0.66, 0.84, 0.15), noxlabels=True),
            add_map_axes(fig, (0.14, 0.48, 0.84, 0.15), noxlabels=True),
            add_map_axes(fig, (0.14, 0.30, 0.84, 0.15), noxlabels=True),
            add_map_axes(fig, (0.14, 0.12, 0.84, 0.15))
        ]
        cx_esa = fig.add_axes((0.14, 0.06, 0.84, 0.015))
    # Column of 4 ESA maps
    elif args.layout == "maps4":
        fig = plt.figure(figsize=(4, 6.5))
        ax_maps = [
            add_map_axes(fig, (0.14, 0.800, 0.84, 0.16), noxlabels=True),
            add_map_axes(fig, (0.14, 0.585, 0.84, 0.16), noxlabels=True),
            add_map_axes(fig, (0.14, 0.370, 0.84, 0.16), noxlabels=True),
            add_map_axes(fig, (0.14, 0.155, 0.84, 0.16))
        ]
        cx_esa = fig.add_axes((0.14, 0.07, 0.84, 0.02))
    # Column of 3 ESA maps
    elif args.layout == "maps3":
        fig = plt.figure(figsize=(4, 5.12))
        ax_maps = [
            add_map_axes(fig, (0.14, 0.75, 0.84, 0.23), noxlabels=True),
            add_map_axes(fig, (0.14, 0.47, 0.84, 0.23), noxlabels=True),
            add_map_axes(fig, (0.14, 0.19, 0.84, 0.23))
        ]
        cx_esa = fig.add_axes((0.14, 0.095, 0.84, 0.023))
    # Hovmöller plot only
    elif args.layout == "hovmoeller":
        fig = plt.figure(figsize=(4, 6.5))
        _, ax_hov = add_hov_axes(
            fig, (0.03, 0.16, 0.79, 0.79),
            latlim=(hovmoeller_extent.stop, hovmoeller_extent.start), # call with S first
            relref=args.target.valid
        )
        cx_esa = fig.add_axes((0.03, 0.08, 0.79, 0.020))
        # Unless forced, set panel letters pool to empty list to suppress "a)"
        # label, as there is only one panel
        if args.panel_letters is None:
            letters = []
    # Scatter plot only
    elif args.layout == "scatter":
        fig = plt.figure(figsize=(4, 4))
        ax_src = None # no scatter source plume panel
        ax_sct = fig.add_axes((0.15, 0.12, 0.80, 0.80))
        if args.panel_letters is None:
            letters = []
    # Cluster-mean maps at scatter valid time
    elif args.layout == "cluster":
        fig = plt.figure(figsize=(4, 4.5))
        ax_top = add_map_axes(fig, (0.12, 0.59, 0.86, 0.30), lonlim=(-135, 75), noxlabels=True)
        ax_bot = add_map_axes(fig, (0.12, 0.24, 0.86, 0.30), lonlim=(-135, 75))
        cx_src = fig.add_axes((0.12, 0.11, 0.86, 0.028))
        tx_src = fig.text(0.12, 0.99, "", ha="left", va="top", fontsize="x-large")
    # Scatter plot with cluster-mean maps at scatter valid time
    elif args.layout == "cluster+scatter":
        fig = plt.figure(figsize=(9, 5))
        ax_src = fig.add_axes((0.68, 0.70, 0.29, 0.24))
        ax_sct = fig.add_axes((0.68, 0.10, 0.29, 0.46))
        ax_top = add_map_axes(fig, (0.06, 0.58, 0.52, 0.30), noxlabels=True)
        ax_bot = add_map_axes(fig, (0.06, 0.21, 0.52, 0.30))
        cx_src = fig.add_axes((0.06, 0.10, 0.52, 0.023))
        tx_src = fig.text(0.06, 0.98, "", ha="left", va="top", fontsize="x-large")
    # Scatter plot with corresponding map on top
    elif args.layout == "map+scatter":
        fig = plt.figure(figsize=(4, 6.3))
        ax_maps = [
            add_map_axes(fig, (0.15, 0.80, 0.80, 0.20))
        ]
        cx_esa = fig.add_axes((0.15, 0.72, 0.80, 0.02))
        ax_src = None # no scatter source plume panel
        ax_sct = fig.add_axes((0.15, 0.08, 0.80, 0.51))
        # The only map shown
        delta_hours = [int((args.scatter.valid - args.target.valid) / dt.timedelta(hours=1))]
    # Scatter plot of F1+F3 and full zonal flux plus linear/parabolic
    # regressions similar to Nakamura and Huang (2018)'s figure 
    elif args.layout == "nh18fig4":
        fig = plt.figure(figsize=(4, 4))
        ax_sct = fig.add_axes((0.18, 0.12, 0.80, 0.80))
        ax_map = add_map_axes(fig, (0.20, 0.70, 0.35, 0.20), lonlim=(-45, 25), latlim=(22, 68), noxlabels=True)
        if args.panel_letters is None:
            letters = []
    # Reanalysis LWA/PV map at target evaluation time and Hovmöllers of LWA and
    # flux similar to Nakamura and Huang (2018)'s figure 5, panels A and C
    elif args.layout == "reanalysis+nh18fig5":
        fig = plt.figure(figsize=(9, 3.5))
        ax_ana = add_map_axes(fig, (0.06, 0.51, 0.48, 0.41), lonlim=(-135, 75))
        _, ax_lwa = add_hov_axes(
            fig, (0.58, 0.08, 0.16, 0.84), lonlim=(-105, 75),
            latlim=(hovmoeller_extent.stop, hovmoeller_extent.start)
        )
        _, ax_flx = add_hov_axes(
            fig, (0.76, 0.08, 0.16, 0.84), lonlim=(-105, 75),
            latlim=(hovmoeller_extent.stop, hovmoeller_extent.start)
        )
        cx_lwa = fig.add_axes((0.06, 0.34, 0.48, 0.03)) # for the LWA colorbar
        cx_flx = fig.add_axes((0.06, 0.14, 0.48, 0.03)) # for the flux colorbar
    # Unknown layout
    else:
        raise NotImplementedError("Layout '{}' is not available".format(args.layout))

    # Remove regions from dataset that are never accessed/visible. This reduces
    # the size of pdf outputs, in which data is included on maps etc. even if
    # it is then clipped. No map shows regions north of 85°N or south of 15°N
    # and no map shows data between 75°E and 165°W (most of Asia and the
    # Pacific).
    def restrict_data_region(d):
        return d.sel(latitude=slice(86, 14), longitude=slice(-164, 76))


    # Read the given datasets containing ensemble data
    datasets = []
    for f in args.data:
        dataset = xr.open_dataset(f)
        # Apply offset so ensemble member numbering is unique
        offset = sum(d.number.size for d in datasets)
        datasets.append(dataset.assign_coords(number=dataset.number.values + offset))
    # Merge all datasets into one
    data = xr.concat(datasets, dim="number", join="inner")
    # Restrict dataset to specified period and reduce pdf output size
    data = data.sel(time=slice(args.start, args.end))
    data = restrict_data_region(data)
    # Resample to the desired frequency
    if args.resample is not None:
        data = data.resample({ "time": args.resample }).mean(keep_attrs=True)
    # Obtain the attributes of the fundamental fields
    pv_attrs  = data.pv.attrs
    lwa_attrs = data.lwa.attrs
    # Basic data integrity checks
    assert "name" in pv_attrs
    assert "unit" in pv_attrs
    assert "level" in pv_attrs
    assert "name" in lwa_attrs
    assert "unit" in lwa_attrs

    # Choose/construct the sensitivity source field
    def _get_source(ds, name):
        if name in (set(ds.variables) - set(ds.coords)):
            return ds[name]
        elif "+" in name:
            terms = name.split("+")
            # Units of all fields must be the same for sum to make sense
            assert len(set(ds[term].attrs["unit"] for term in terms)) == 1, "incompatible units of source fields in sum"
            # New DataArray as sum of fields
            out = reduce(operator.add, (ds[term] for term in terms))
            out.attrs["name"] = "+".join(ds[term].attrs["name"] for term in terms)
            out.attrs["unit"] = ds[terms[0]].attrs["unit"]
            return out
        else:
            raise ValueError("invalid source field '{}'".format(args.source))
    # Obtain source field from the ensemble
    source_ens = _get_source(data, args.source)

    # Initialize the target metric
    if args.target_metric == "average":
        target_metric = metrics.BoxAverage(args.target)
    elif args.target_metric.startswith("area>"):
        _target_threshold = float(args.target_metric[5:])
        target_metric = metrics.BoxThresholdArea(args.target, _target_threshold)
    else:
        raise ValueError("Invalid target metric '{}'".format(args.target_metric))

    # Evaluate target metric for ensemble data
    target_ens = target_metric(data.lwa)
    target_ens_values = target_ens.sel(time=args.target.valid).values
    # Compute rank-sorted indices of members for cluster selection
    target_ens_ranks = np.argsort(target_ens_values)
    
    # Member indices for clusters
    bot_cluster = target_ens_ranks[              : args.ncluster]
    not_cluster = target_ens_ranks[ args.ncluster:-args.ncluster]
    top_cluster = target_ens_ranks[-args.ncluster:              ]

    # Same for reanalysis data if given
    if args.reanalysis is not None:
        # Read given reanalysis dataset, restrict to specified period, resample
        ana = xr.open_dataset(args.reanalysis)
        ana = ana.sel(time=slice(args.start, args.end))
        ana = restrict_data_region(ana)
        if args.resample is not None:
            ana = ana.resample({ "time": args.resample }).mean(keep_attrs=True)
        # Extract target metric from reanalysis (restrict to time period
        # covered by ensemble)
        target_ana = target_metric(ana.lwa).sel(time=slice(min(data.time.values), max(data.time.values)))
        target_ana_values = target_ana.sel(time=args.target.valid).values
        # Must have same attributes as ensemble data
        assert pv_attrs["name"] == ana.pv.attrs["name"]
        assert pv_attrs["level"] == ana.pv.attrs["level"]
        assert pv_attrs["unit"] == ana.pv.attrs["unit"]
        assert lwa_attrs["name"] == ana.lwa.attrs["name"]
        assert lwa_attrs["unit"] == ana.lwa.attrs["unit"]
        # Obtain source field from the reanalysis
        source_ana = _get_source(ana, args.source)


    # Framework-specific Configuration/Plotters #

    ISQG = (data.attrs["framework"] == "QG")
    lwa_levels = np.arange(0, 1000, (20 if ISQG else 50))
    flux_levels = [
        -2000, -1800, -1600, -1400, -1200,
        -1000,  -800,  -600,  -400,  -200,
          200,   400,   600,   800,  1000,
         1200,  1400,  1600,  1800,  2000
    ]

    def contour_pv(ax, y, x, field, linewidths=1.2):
        return plotting.contour(ax, y, x, field, **{
            "levels": [1.0e-4, 2.0e-4] if ISQG else [2.0e-6, 4.0e-6],
            "linewidths": linewidths
        })

    def contourf_lwa(ax, y, x, field):
        return plotting.contourf_max(ax, y, x, field, **{
            "levels": lwa_levels[:len(plotting.lwa_colors)],
            "colors": plotting.lwa_colors
        })

    def contourf_flux(ax, y, x, field):
        k = (len(flux_levels) - 1 - len(plotting.flux_colors)) // 2
        return plotting.contourf_sym(ax, y, x, field, **{
            "levels": flux_levels[k+1:k+len(plotting.flux_colors)],
            "colors": plotting.flux_colors
        })

    def contourf_source(ax, y, x, field):
        if args.source == "lwa":
            return contourf_lwa(ax, y, x, field)
        elif "f1" in args.source or "f2" in args.source or "f3" in args.source:
            return contourf_flux(ax, y, x, field)
        else:
            raise ValueError("No contourf preset for source field '{}'".format(args.source))

    def contours_esa(ax, y, x, src, tgt, cluster=None):
        # Source-mean as contours, optionally only for a given cluster
        field = (src if cluster is None else src[cluster,...]).mean(axis=0)
        ct = plotting.contour(ax, y, x, field, **{
            "levels": lwa_levels if args.source == "lwa" else flux_levels,
            "linewidths": 0.8
        })
        plt.clabel(ct, ct.levels[::2], fontsize="x-small")
        # Sensitivity field as filled contours
        return plotting.contourf_corr(ax, y, x, src, tgt, **{
            "nsamples": args.nsamples,
            "sample_size": args.sample_size,
            "significance_level": args.significance_level
        })

    # Use colors from clustering for scatter plots
    sct_color = np.array(["#333333"] * data.number.size)
    sct_color[top_cluster] = esa_red
    sct_color[bot_cluster] = esa_blue

    # Reanalysis Panel #

    if "reanalysis" in args.layout:
        assert args.reanalysis is not None, "No reanalysis data given for reanalysis panel"
        # Plot LWA and PV in the reanalysis panel
        contourf_lwa(
            ax_ana,
            ana.latitude.values,
            ana.longitude.values,
            ana.lwa.sel(time=args.target.valid).values
        )
        ct = contour_pv(
            ax_ana,
            ana.latitude.values,
            ana.longitude.values,
            ana.pv.sel(time=args.target.valid).values,
            linewidths=1.6
        )
        # Configure reanalysis panel
        set_panel_title(ax_ana, "{pv_name} ({pv_level}), {lwa_name}".format(
            lwa_name=texify(lwa_attrs["name"]),
            pv_name=texify(pv_attrs["name"]),
            pv_level=texify(pv_attrs["level"])
        ))
        ax_ana.set_title(format_date_full(args.target.valid), loc="right")
        plot_box(ax_ana, args.target, linewidth=2)


    # Plume and Rank Cluster Panels #

    if "plume" in args.layout:
        # Vertical line marking the target valid time
        ax_plu.axvline(args.target.valid, color="#000000", linewidth=1.5, linestyle="dashed")
        # Plot timeseries of rank cluster and other members
        _legend = plotting.plume(
            ax_plu,
            target_ens.time.values,
            target_ens.values,
            clusters=(bot_cluster, not_cluster, top_cluster),
            reanalysis=(None if args.reanalysis is None else target_ana.values),
            ncluster=args.ncluster
        )
        ax_plu.legend(*_legend, loc="upper left", fontsize="small")
        # Only use 00 times as time axis labels
        time_ticks = [t.values for t in target_ens.time if t.dt.hour == 0]
        ax_plu.set_xticks(time_ticks[::2])
        ax_plu.set_xticks(time_ticks, minor=True)
        ax_plu.xaxis.set_major_formatter(date_formatter)
        ax_plu.set_ylim(mult_limit(target_ens.values, mult=5, tol=5, minval=0))
        ax_plu.set_xlim((min(target_ens.time.values), max(target_ens.time.values)))
        set_grid(ax_plu)
        ax_plu.set_ylabel(target_metric.label(lwa_attrs))
        set_panel_title(ax_plu, "Evolution of Target Metric")

        # Cluster-mean LWA and PV for top members
        cf_lwa = contourf_lwa(
            ax_top,
            data.latitude.values,
            data.longitude.values,
            data.lwa.sel(time=args.target.valid).values[top_cluster,:,:].mean(axis=0)
        )
        contour_pv(
            ax_top,
            data.latitude.values,
            data.longitude.values,
            data.pv.sel(time=args.target.valid).values[top_cluster,:,:].mean(axis=0)
        )
        # Configure top rank plot
        _title = "Top Cluster {} ({}), {}".format(
            texify(pv_attrs["name"]), texify(pv_attrs["level"]), texify(lwa_attrs["name"])
        )
        set_panel_title(ax_top, _title, fontsize="medium", color=esa_red)
        ax_top.spines["geo"].set_color(esa_red)
        ax_top.tick_params(axis="both", color=esa_red)
        plot_box(ax_top, args.target, linewidth=2)

        # Extract the LWA colorbar from this panel
        cb = plt.colorbar(cf_lwa, cax=cx_lwa)
        cb.set_label("{name} [{unit}]".format(**texify(lwa_attrs)))

        # Cluster-mean LWA and PV for bottom members
        contourf_lwa(
            ax_bot,
            data.latitude.values,
            data.longitude.values,
            data.lwa.sel(time=args.target.valid).values[bot_cluster,:,:].mean(axis=0)
        )
        contour_pv(
            ax_bot,
            data.latitude.values,
            data.longitude.values,
            data.pv.sel(time=args.target.valid).values[bot_cluster,:,:].mean(axis=0)
        )
        # Configure bottom rank plot
        _title = "Bottom Cluster {} ({}), {}".format(
            texify(pv_attrs["name"]), texify(pv_attrs["level"]), texify(lwa_attrs["name"])
        )
        set_panel_title(ax_bot, _title, fontsize="medium", color=esa_blue)
        ax_bot.spines["geo"].set_color(esa_blue)
        ax_bot.tick_params(axis="both", color=esa_blue)
        plot_box(ax_bot, args.target, linewidth=2)


    # ESA Map Panels #

    if "map" in args.layout:
        # Time deltas for maps were specified above
        for ax_map, delta_hour in zip(ax_maps, delta_hours):
            valid = args.target.valid + dt.timedelta(hours=delta_hour)
            # Plot correlation map and LWA contours of top cluster
            cf_esa = contours_esa(
                ax_map,
                data.latitude.values,
                data.longitude.values,
                source_ens.sel(time=valid).values,
                target_ens_values,
                cluster=top_cluster
            )
            # Show the target box in the delta = 0 panel, but only if the
            # target field is the same as the source field
            if delta_hour == 0 and args.source == "lwa":
                plot_box(ax_map, args.target, linewidth=2)
            # Show the scatter box if specified and the valid times match
            if args.scatter is not None and valid == args.scatter.valid:
                plot_box(ax_map, args.scatter, linewidth=2, linestyle="solid")
            # Configure map panel
            set_panel_title(ax_map, "ESA: {name} [{unit}]".format(**texify(source_ens.attrs)))
            ax_map.set_title("{:+} h".format(delta_hour), loc="right")
        # Correlation colorbar
        cb = plt.colorbar(cf_esa, cax=cx_esa, orientation="horizontal", spacing="proportional")
        cb.set_label("Sensitivity to {name} (Correlation)".format(**texify(source_ens.attrs)))


    # Hovmöller Panel #

    if "hovmoeller" in args.layout:
        # Reduce data along meridional dimension
        src_slice = source_ens.sel(latitude=hovmoeller_extent)
        weights = np.cos(np.deg2rad(src_slice.latitude))
        # Plot correlation map and LWA contours of top cluster
        cf_esa = contours_esa(
            ax_hov,
            data.time.values,
            data.longitude.values,
            src_slice.weighted(weights).mean(["latitude"]).values.swapaxes(0, 1),
            target_ens_values,
            cluster=top_cluster
        )
        # Correlation colorbar
        cb = plt.colorbar(cf_esa, cax=cx_esa, orientation="horizontal", spacing="proportional")
        cb.set_label("Sensitivity to {name} (Correlation)".format(**texify(source_ens.attrs)))
        # Line to indicate target valid time
        ax_hov.axhline(args.target.valid, linestyle="dashed", linewidth=1.5, color="#000000")
        # Configure Hovmöller plot
        set_panel_title(ax_hov, "ESA: {name} [{unit}]".format(**texify(source_ens.attrs)))
        ax_hov.set_title("{:%Y}".format(args.target.valid), loc="right")


    # Cluster Map Panels #

    if "cluster" in args.layout:
        # If scatter parameters are set, use scatter valid time, else target
        _valid = args.target.valid if args.scatter is None else args.scatter.valid
        # Source field information in main caption
        tx_src.set_text(format_date_full(_valid))
        # Plot means of PV and source field for both clusters
        for _ax, _cluster in [(ax_top, top_cluster), (ax_bot, bot_cluster)]:
            cf_src = contourf_source(
                _ax,
                data.latitude.values,
                data.longitude.values,
                source_ens.sel(time=_valid).values[_cluster,:,:].mean(axis=0)
            )
            contour_pv(
                _ax,
                data.latitude.values,
                data.longitude.values,
                data.pv.sel(time=_valid).values[_cluster,:,:].mean(axis=0)
            )
            if args.scatter is not None:
                plot_box(_ax, args.scatter, linewidth=2, linestyle="solid")
        # Colorbar for source field
        cb = plt.colorbar(cf_src, cax=cx_src, orientation="horizontal", spacing="proportional")
        cb.set_label("{name} [{unit}]".format(**texify(source_ens.attrs)))
        # Titles and coloring of cluster panels
        _title = "Top Cluster {} ({}), {}".format(
            texify(pv_attrs["name"]), texify(pv_attrs["level"]), texify(source_ens.attrs["name"])
        )
        set_panel_title(ax_top, _title, color=esa_red)
        ax_top.spines["geo"].set_color(esa_red)
        _title = "Bottom Cluster {} ({}), {}".format(
            texify(pv_attrs["name"]), texify(pv_attrs["level"]), texify(source_ens.attrs["name"])
        )
        set_panel_title(ax_bot, _title, color=esa_blue)
        ax_bot.spines["geo"].set_color(esa_blue)


    # Scatter Panel #

    if "scatter" in args.layout:
        assert args.scatter is not None, "No scatter extraction valid time and box given"
        _scatter_metric = metrics.BoxAverage(args.scatter)
        # Show plume of source metric for 4 days around scatter source valid
        # time if a corresponding panel exists
        if ax_src is not None:
            _period = { "time": slice(
                    args.scatter.valid - dt.timedelta(hours=48),
                    args.scatter.valid + dt.timedelta(hours=48)
            ) }
            source_ens_plume = _scatter_metric(source_ens.sel(_period))
            # Mark valid time of scatter source
            ax_src.axvline(args.scatter.valid, color="#000000", linewidth=1.5, linestyle="dashed")
            # Plot timeseries of rank cluster and other members, add the
            # reanalysis timeseries if data is available
            _legend = plotting.plume(
                ax_src,
                source_ens_plume.time.values,
                source_ens_plume.values,
                clusters=(bot_cluster, not_cluster, top_cluster),
                reanalysis=(None if args.reanalysis is None else _scatter_metric(source_ana.sel(_period)).values),
                linewidth=0.8,
                ncluster=args.ncluster
            )
            ax_src.legend(*_legend, loc="upper left", fontsize="x-small", handlelength=1.0)
            # Only use 00 times as time axis labels
            ax_src.set_xticks([t.values for t in source_ens_plume.time if t.dt.hour == 0])
            ax_src.set_xticks([t.values for t in source_ens_plume.time if t.dt.hour == 12], minor=True)
            ax_src.xaxis.set_major_formatter(date_formatter)
            ax_src.set_xlim((min(source_ens_plume.time.values), max(source_ens_plume.time.values)))
            ax_src.set_ylim(mult_limit(source_ens_plume.values[:-2,:], 50))
            set_panel_title(ax_src, "Evolution of Agg. " + texify(source_ens.attrs["name"]))
            ax_src.set_ylabel("Agg. {source} [{unit}]".format(
                source=texify(source_ens.attrs["name"]),
                unit=texify(source_ens.attrs["unit"])
            ))
        # Evaluate at scatter source valid time
        source_ens_values = _scatter_metric(source_ens.sel(time=args.scatter.valid)).values
        # Markers for ENS
        ax_sct.scatter(
            source_ens_values,
            target_ens_values,
            s=25,
            c=sct_color,
            marker="o",
            alpha=0.80,
            linewidth=0,
            label="ENS member"
        )
        # Add marker for reanalysis if available
        if args.reanalysis is not None:
            ax_sct.scatter(
                _scatter_metric(source_ana.sel(time=args.scatter.valid)).values,
                target_ana_values,
                s=125,
                c="#000000",
                marker="*",
                label="Reanalysis"
            )
        # Linear fit to ensemble
        reg = linregress(source_ens_values, target_ens_values)
        # Idealized reference curve if data given
        if args.scatter_idealized is not None:
            with open(args.scatter_idealized, "r") as f:
                idealized = { k: np.asarray(v) for k, v in json.load(f).items() }
            _a = idealized["z"] # jammed 
            _b = ~idealized["z"] | np.roll(~idealized["z"], 1) # non-jammed
            ax_sct.plot(idealized["x"][_a], idealized["y"][_a], color="#000000", linewidth=1.8, linestyle="solid", label="Idealized")
            ax_sct.plot(idealized["x"][_b], idealized["y"][_b], color="#000000", linewidth=1.8, linestyle=(0, (2, 1)))
        else:
            xmin = min(source_ens_values) - 20
            xmax = max(source_ens_values) + 20
            ax_sct.plot(
                [xmin, xmax],
                [reg.intercept + reg.slope * xmin, reg.intercept + reg.slope * xmax],
                color="#000000",
                linewidth=1.
            )
        # Configure scatter panel
        ax_sct.legend(loc="upper left", fontsize="x-small", markerscale=0.7, handlelength=1.0)
        ax_sct.set_ylim(mult_limit(target_ens_values, mult=5, minval=0))
        ax_sct.set_ylabel("Target ({}) [{}]".format(
            format_date_full(args.target.valid),
            texify(lwa_attrs["unit"])
        ))
        ax_sct.set_xlim(mult_limit(source_ens_values, mult=100, tol=50))
        ax_sct.set_xlabel("Aggregated {source} ({dt:+} h) [{unit}]".format(
            source=texify(source_ens.attrs["name"]),
            dt=int((args.scatter.valid - args.target.valid) / dt.timedelta(hours=1)),
            unit=texify(source_ens.attrs["unit"])
        ))
        set_grid(ax_sct)
        set_panel_title(ax_sct, "Sensitivity Scatter")
        ax_sct.set_title("$r = {:.2f}$".format(reg.rvalue), loc="right", fontsize="medium")


    # Nakamura and Huang (2018) figure 4 #

    if "nh18fig4" in args.layout:
        # The variables shown in this plot are fixed and specially aggregated
        # to match the processing of NH18 (in particular a 4 day temporal
        # average around the scatter valid time is computed).
        assert args.scatter is not None, "No scatter extraction valid time and box given"
        # Extract a 4-day period, left inclusive, right exclusive
        _tslc = slice(
            args.scatter.valid - dt.timedelta(hours=48),
            args.scatter.valid + dt.timedelta(hours=47)
        )
        # Instead of extracting a single gridpoint, use a box
        _metric = metrics.BoxAverage(args.scatter)
        # Extract the required fields from the ENS: F1+F3, total F, lwa
        _ens4d = data.sel(time=_tslc)
        _f1f3_ens = _metric(_get_source(_ens4d, "f1+f3").mean(dim="time")).values
        _ftot_ens = _metric(_get_source(_ens4d, "f1+f2+f3").mean(dim="time")).values
        _lwa_ens  = _metric(_ens4d.lwa.mean(dim="time")).values
        # And the same from the reanalysis
        _ana4d = ana.sel(time=_tslc)
        _f1f3_ana = _metric(_get_source(_ana4d, "f1+f3").mean(dim="time")).values
        _ftot_ana = _metric(_get_source(_ana4d, "f1+f2+f3").mean(dim="time")).values
        _lwa_ana  = _metric(_ana4d.lwa.mean(dim="time")).values
        # Line plots of NH18 fits are omitted in the public code
        # Round markers for the ensemble members
        ax_sct.scatter(_lwa_ens, _f1f3_ens, s=25, marker="o", alpha=0.8, linewidth=0, c="#F90", label=texify("F1+F3"))
        ax_sct.scatter(_lwa_ens, _ftot_ens, s=25, marker="o", alpha=0.8, linewidth=0, c=sct_color, label=texify("F1+F2+F3"))
        # A star each for the reanalysis values
        ax_sct.scatter(_lwa_ana, _f1f3_ana, s=125, marker="*", c="#C60")
        ax_sct.scatter(_lwa_ana, _ftot_ana, s=125, marker="*", c="#000")
        # Style panel
        ax_sct.set_xlim((25, 95))
        ax_sct.set_xlabel("{name} [{unit}]".format(**texify(lwa_attrs)))
        ax_sct.set_ylim((-500, 2000))
        ax_sct.set_ylabel("Flux [{unit}]".format(**texify(data.f1.attrs)))
        ax_sct.legend(loc="lower center", ncol=2, handlelength=1.0)
        set_grid(ax_sct)
        set_panel_title(ax_sct, "{} – {} {:%Y}".format(
            date_formatter.format_data(_tslc.start),
            date_formatter.format_data(_tslc.stop),
            _tslc.stop
        ))
        ax_sct.set_title("ENS, ERA5", loc="right")
        # Add the averaging region in the map panel
        plot_box(ax_map, args.scatter, linewidth=2, linestyle="solid", fill=True, facecolor="#F7F7F7")
        ax_map.yaxis.tick_right()
        ax_map.set_yticks([])
        ax_map.set_xticks([])
        # Print the correlation value for F1+F3 vs. <A>cos(ϕ)
        print("corr[<A>cos(ϕ), F1+F3] = {:.2f}".format(linregress(_lwa_ens, _f1f3_ens).rvalue))


    # LWA and flux Hovmöller panels of reanalysis #

    if "nh18fig5" in args.layout:
        assert args.reanalysis is not None
        # Displayed variables are fixed and do not depend on source
        _lwa = ana.lwa.sel(latitude=hovmoeller_extent)
        _flx = _get_source(ana, "f1+f2+f3").sel(latitude=hovmoeller_extent)
        weights = np.cos(np.deg2rad(_lwa.latitude))
        # Hovmöller fields of LWA and Flux
        _lwa_args = (
            _lwa.time.values,
            _lwa.longitude.values,
            _lwa.weighted(weights).mean(["latitude"]).values
        )
        _flx_args = (
            _flx.time.values,
            _flx.longitude.values,
            _flx.weighted(weights).mean(["latitude"]).values
        )
        # 60 m/s contour of LWA as in NH18
        plt.clabel(plotting.contour(ax_lwa, *_lwa_args, levels=[60]), fontsize="x-small")
        # LWA field with usual colorbar
        cf_lwa = contourf_lwa(ax_lwa, *_lwa_args)
        cb_lwa = plt.colorbar(cf_lwa, cax=cx_lwa, orientation="horizontal", spacing="proportional")
        cb_lwa.set_label("{name} [{unit}]".format(**texify(_lwa.attrs)))
        # 0 and 600 m²/s² contours of flux as in NH18
        plt.clabel(plotting.contour(ax_flx, *_flx_args, levels=[0, 600], colors=["#999", "#000"]), fontsize="x-small")
        # Flux field with usual colorbar
        cf_flx = contourf_flux(ax_flx, *_flx_args)
        cb_flx = plt.colorbar(cf_flx, cax=cx_flx, orientation="horizontal", spacing="proportional")
        cb_flx.set_label("{name} [{unit}]".format(**texify(_flx.attrs)))
        # Style panels
        ax_lwa.set_yticklabels([])
        for _ax in [ax_lwa, ax_flx]:
            _ax.axhline(args.target.valid, color="k", linewidth=1.5, linestyle="dashed")
        set_panel_title(ax_lwa, "{name}".format(**texify(_lwa.attrs)))
        set_panel_title(ax_flx, "{name}".format(**texify(_flx.attrs)))


    # Output #

    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output, dpi=200)

