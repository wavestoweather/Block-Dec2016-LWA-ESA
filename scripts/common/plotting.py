import numbers
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from cartopy.mpl.geoaxes import GeoAxes

from .esa import esa_corr_with_confidence


# Custom colors and colormaps
lwa_colors  = ["#FFFFFF", "#C9F0C4", "#84DEBF", "#66B2DA", "#7972D8", "#903C9A", "#7B2140", "#3b1705"]
flux_colors = ["#c51b7d", "#de77ae", "#f1b6da", "#fde0ef", "#f7f7f7", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221"]
esa_colors  = ["#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"]
esa_levels = [-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8] # use with extend="both"
esa_blue = esa_colors[ 0] # bottom rank cluster
esa_red  = esa_colors[-1] # top rank cluster


# Instead of using grey for the grid, use semi-transparent black lines
grid_rgb = "#000000"
grid_alpha = 0.2
grid_rgba = "{}{:0>2}".format(grid_rgb, hex(int(255 * grid_alpha))[2:])

def set_grid(ax, **kwargs):
    ax.grid(color=grid_rgb, alpha=grid_alpha)


# Coordinate system configuration and labelling
lat_ticks  = [20, 40, 60, 80]
lat_labels = ["20°N", "40°N", "60°N", "80°N"]
lon_ticks  = [-150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150]
lon_labels = ["150°W", "120°W", "90°W", "60°W", "30°W", "0°", "30°E", "60°E", "90°E", "120°E", "150°E"]
lon_labels_sparse = ["", "120°W", "", "60°W", "", "0°", "", "60°E", "", "120°E", ""]

# Create PlateCarree maps
projection = ccrs.PlateCarree()
# Data is given on lat-lon grid
transform  = ccrs.PlateCarree()

# Date and time labels for axes, show only the day as a number and the month as
# abbreviated text
def date_formatter_func(value, pos, relref=None):
    if isinstance(value, numbers.Number):
        value = mdates.num2date(value).replace(tzinfo=None)
    # Cannot rely on %b as it depends on locale, use own list of month abbrevs
    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    if relref is None:
        return "{} {}".format(value.day, months[value.month-1])
    else:
        return "{day} {month}\n({dt:+} h)".format(
            day=value.day,
            month=months[value.month-1],
            dt=int((value - relref) / dt.timedelta(hours=1))
        )
date_formatter = mticker.FuncFormatter(date_formatter_func)

# An associated format function that also includes the year and hour
def format_date_full(x):
    return "{} {:%Y %H} UTC".format(date_formatter.format_data(x), x)

# Axes limits as a multiple of some value
def mult_limit(values, mult=10, tol=0, minval=None, maxval=None):
    # Find next multiple below minium
    ymin = int(min(values) / mult) * mult
    while min(values) < ymin:
        ymin -= mult
    # Find mext multiple above maximum
    ymax = int(max(values) / mult) * mult
    while ymax < max(values):
        ymax += mult
    # Postprocessing: add tolerance
    if tol < 0:
        tol = -tol * (ymax - ymin)
    ymin = ymin - tol
    ymax = ymax + tol
    # Postprocessing: restrict to given bounds
    if minval is not None and ymin < minval:
        ymin = minval
    if maxval is not None and ymax > maxval:
        ymax = maxval
    return ymin, ymax

# Add a pre-configured GeoAxes to a Figure
def add_map_axes(fig, rect, latlim=(15, 85), lonlim=(-165, 75), noxlabels=False, **kwargs):
    ax = fig.add_axes(rect, projection=projection, anchor="SW", **kwargs)
    ax.coastlines(linewidth=0.5, color="#555555")
    grid = ax.gridlines(
        xlocs=lon_ticks,
        ylocs=lat_ticks,
        draw_labels=False,
        linewidth=0.5,
        color=grid_rgba
    )
    ax.set_xticks(lon_ticks)
    if noxlabels:
        xticklabels = []
    elif lonlim[1] - lonlim[0] < 180:
        xticklabels = lon_labels
    else:
        xticklabels = lon_labels_sparse
    ax.set_xticklabels(xticklabels)
    ax.set_yticks(lat_ticks)
    ax.set_yticklabels(lat_labels)
    # Set limits at the end otherwise ticks overwrite
    ax.set_xlim(lonlim)
    ax.set_ylim(latlim)
    return ax

# Add a pre-configured longitude-time Hovmöller Axes to a Figure
def add_hov_axes(fig, rect, latlim, lonlim=(-165, 75), relref=None, **kwargs):
    # Use a map of the reduction-region as the x-labels
    ax_map = add_map_axes(fig, rect, lonlim=lonlim, latlim=latlim, **kwargs)
    ax_map.set_yticks([])
    # The cartopy GeoAxes will generally not fill the entire height of the box
    # specified for ax_map, extract the real height and adapt the Hovmöller
    # panel to be flush with the map (which is SW-anchored)
    ax_map_top = ax_map.get_position().y1
    ax_hov = fig.add_axes((rect[0], ax_map_top, rect[2], rect[1] + rect[3] - ax_map_top))
    set_grid(ax_hov) 
    ax_hov.set_xticks(lon_ticks)
    ax_hov.xaxis.set_tick_params(bottom=False, labelbottom=False)
    ax_hov.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda value, pos: date_formatter_func(value, pos, relref=relref))
    )
    ax_hov.yaxis.tick_right()
    ax_hov.set_xlim(lonlim)
    return ax_map, ax_hov


def get_transform(ax):
    return transform if isinstance(ax, GeoAxes) else ax.transData

# Show a (target) box with an outline on a map
def plot_box(ax, box, **kwargs):
    props = {
        "fill": False,
        "edgecolor": "black",
        "linestyle": "dashed",
        "transform": get_transform(ax)
    }
    props.update(kwargs)
    width  = box.E - box.W
    height = box.N - box.S
    rect = plt.Rectangle([box.W, box.S], width, height, **props)
    ax.add_patch(rect)


# Plot presets

def contour(ax, y, x, field, **kwargs):
    props = {
        "linewidths": 1,
        "colors": "#000000",
        "transform": get_transform(ax)
    }
    props.update(kwargs)
    return ax.contour(x, y, field, **props)

def contourf(ax, y, x, field, **kwargs):
    props = {
        "transform": get_transform(ax)
    }
    props.update(kwargs)
    return ax.contourf(x, y, field, **props)

def contourf_max(ax, y, x, field, **kwargs):
    props = {
        "colors": "cubehelix_r",
        "extend": "max"
    }
    props.update(kwargs)
    return contourf(ax, y, x, field, **props)

def contourf_sym(ax, y, x, field, **kwargs):
    props = {
        "colors": "RdBu_r",
        "extend": "both"
    }
    props.update(kwargs)
    return contourf(ax, y, x, field, **props)

# Sensitivity field filled contours and optional non-significant regions
def contourf_corr(ax, y, x, src, tgt, nsamples=0, sample_size=None, significance_level=0.1, **kwargs):
    correlation, confidence = esa_corr_with_confidence(src, tgt, nsamples=nsamples, sample_size=sample_size)
    # Correlation map
    props = {
        "levels": esa_levels,
        "colors": esa_colors,
        "extend": "both",
        "transform": get_transform(ax)
    }
    props.update(kwargs)
    cf = ax.contourf(x, y, correlation, **props)
    # Stippling of significant regions if confidence intervals are given
    if confidence is not None:
        # Matplotlib's build-in hatching is quite limited, try to do it better
        # with scatter. Need to calculate good spacing for the points.
        figw, figh = ax.get_figure().get_size_inches()
        _, _, axw, axh = ax.get_position().bounds
        ratio = ax.get_data_ratio() * (axw * figw) / (axh * figh)
        # Control the density of points, reuse the data grid (good enough)
        dx = round(x.size / 40) # controls the density of points
        dy = max(1, round(dx * ratio))
        slc = slice(dy // 2, None, dy), slice(dx // 2, None, dx)
        # Give the points a red or blue tint to make them blend better with the
        # underlying filled contours
        is_significant = (confidence <= significance_level)[slc]
        confidence_slc = confidence[slc][is_significant]
        colors = np.full(confidence_slc.shape, "#333333") # grey
        colors[correlation[slc][is_significant] >  0.2] = "#650E19" # dark red
        colors[correlation[slc][is_significant] < -0.2] = "#123960" # dark blue
        # Place points
        xx, yy = np.meshgrid(x, y)
        ax.scatter(xx[slc][is_significant], yy[slc][is_significant], s=1, c=colors, transform=get_transform(ax))
    # Return for colorbar setup
    return cf

