import argparse
import datetime as dt
import re
import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from .common import metrics, plotting


def verify_path(arg):
    assert os.path.exists(arg), "'{}' does not exist".format(arg)
    return arg

def verify_metric(arg):
    assert arg == "average" or arg.startswith("area>")
    return arg

def parse_region(arg):
    N, E, S, W = [float(x) for x in arg.split(",")]
    if N < S:
        N, S = S, N
    if E < W:
        E, W = W, E
    return metrics.TargetBox(N, E, S, W)

def parse_datehour(arg):
    return dt.datetime.strptime(arg, "%Y-%m-%dT%H")

parser = argparse.ArgumentParser()
parser.add_argument("--reanalysis", type=verify_path, default=None)
parser.add_argument("--target-metric", type=verify_metric, default="average")
parser.add_argument("--output", type=str, default=None)
parser.add_argument("--highlight", type=parse_datehour, nargs="*", default=[], metavar="yyyy-mm-ddThh")
parser.add_argument("--percentile", type=int, default=10, help="flyer percentile (default 10)")
parser.add_argument("target", type=parse_region, metavar="N,E,S,W", help="target region")
parser.add_argument("data", type=verify_path, help="input forecast data")

if __name__ == "__main__":
    args = parser.parse_args()

    # Initialize the target metric
    if args.target_metric == "average":
        target_metric = metrics.BoxAverage(args.target)
    elif args.target_metric.startswith("area>"):
        _target_threshold = float(args.target_metric[5:])
        target_metric = metrics.BoxThresholdArea(args.target, _target_threshold)
    else:
        raise ValueError("Invalid target metric '{}'".format(args.target_metric))

    # Process forecast data
    ens = xr.open_dataset(args.data)
    target_ens = target_metric(ens.lwa).to_dataframe().unstack()
    target_valid = pd.to_datetime(ens.valid.values)[0]

    # Process reanalysis data if given
    if args.reanalysis is not None:
        ana = xr.open_dataset(args.reanalysis)
        ana = ana.sel(time=target_valid)
        target_ana = target_metric(ana.lwa)

    # Single-column figure
    fig = plt.figure(figsize=(4, 2.5))

    # Boxplot axes
    ax = fig.add_axes((0.13, 0.15, 0.84, 0.74))
    ax.set_title("Target Metric Forecasts", loc="left")

    # Draw boxplots individually so they can be highlighted
    for i, (init, values) in enumerate(target_ens.iterrows()):
        color = "#000" if init in args.highlight else "#444"
        bp = ax.boxplot(
            values,
            positions=[i],
            whis=(args.percentile, 100-args.percentile),
            patch_artist=True,
            widths=0.7,
            boxprops={ "color": color, "facecolor": color },
            whiskerprops={ "linewidth": 3.5, "color": color },
            showcaps=False,
            medianprops={ "color": "#999" },
            flierprops={ "visible": False }
        )
        # Take positions of (invisible) fliers and re-plot them with
        # appropriate cluster colors
        fx, fy = bp["fliers"][0].get_data()
        half = len(fx) // 2
        ax.plot(fx[:half], fy[:half], linewidth=0, color=plotting.esa_blue, marker="x", markersize=4, zorder=10)
        ax.plot(fx[half:], fy[half:], linewidth=0, color=plotting.esa_red, marker="x", markersize=4, zorder=10)

    # Add marker and horizontal line for reanalysis value if available
    if args.reanalysis is not None:
        assert init == target_valid
        ax.plot([i], [target_ana], marker="*", markersize=10, color="#000", linewidth=1, zorder=99, label="Reanalysis")
        ax.axhline(target_ana, linewidth=1, color="#000", zorder=0)
        ax.legend(loc="center right", fontsize="small")

    # Translate x-ticks to dates and add nice labels
    ticks = np.arange(0, len(target_ens)+1, 4)
    ax.set_xticks(ticks)
    ax.set_xticks(np.arange(0, len(target_ens)+1, 2), minor=True)
    def init_label(d):
        out = plotting.date_formatter.format_data(d)
        lead = int((target_valid - d).total_seconds() / 60 / 60)
        if lead == 0:
            return out + "\n{:%Y}".format(d)
        else:
            return out + "\n(-{} h)".format(lead)
    ax.set_xticklabels(target_ens.index[ticks].map(init_label), fontsize="small")
    ax.set_xlim((-1, len(target_ens)))

    # Use default y-ticks, only add label
    ax.set_ylabel(target_metric.label(ens.lwa.attrs), fontsize="small")

    # Add a small map in the bottom right showing the target region
    mx = plotting.add_map_axes(fig, (0.705, 0.19, 0.35, 0.26), lonlim=(-60, 30), latlim=(20, 80), noxlabels=True)
    mx.set_xticks([])
    mx.set_yticks([])
    plotting.plot_box(mx, args.target, linewidth=2)

    # Show or write to disk
    if args.output is None:
        plt.show()
    else:
        fig.savefig(args.output, dpi=200)

