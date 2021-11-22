import collections
import numpy as np

from .texify import texify


class _BoxMetricBase:

    def _restrict_region(self, da):
        return da.sel(
            latitude=slice(self.box.N, self.box.S),
            longitude=slice(self.box.W, self.box.E)
        )


TargetBox = collections.namedtuple("TargetProperties", ["N", "E", "S", "W"])

class BoxAverage(_BoxMetricBase):

    def __init__(self, box):
        self.box = box

    def __call__(self, da):
        region = self._restrict_region(da)
        # https://xarray.pydata.org/en/stable/examples/area_weighted_temperature.html
        weights = np.cos(np.deg2rad(region.latitude))
        return region.weighted(weights).mean(["latitude", "longitude"])

    def label(self, field_attrs):
        return "Averaged {name} [{unit}]".format(**texify(field_attrs))


class BoxThresholdArea(_BoxMetricBase):

    def __init__(self, box, threshold):
        self.box = box
        self.threshold = threshold

    def __call__(self, da):
        region = self._restrict_region(da)
        weights = np.cos(np.deg2rad(region.latitude))
        # Total area in region for normalization
        total = weights.sum() * region.longitude.size
        # Normalize with area of full box to scale into 0% to 100% interval
        out = (weights * (region >= self.threshold)).sum(["latitude", "longitude"]) / total * 100.
        out.name = "area_ge_{}".format(self.threshold)
        return out

    def label(self, field_attrs):
        return "{name} > {threshold:.0f} {unit} Fraction [%]".format(
            threshold=self.threshold, **texify(field_attrs)
        )

