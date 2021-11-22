import numpy as np
import cffi
from . import _regrid


# CFFI preparation
_ffi = cffi.FFI()

def _double_pointer(arr):
    return _ffi.cast("double *", _ffi.from_buffer(arr))


# Bindings to compiled regridding to half resolution
def half_resolution(lat, lon, *fields):
    """Half the resolution of the given input fields, with an area-weighted
    interpolation method. This function only works for grids with an odd number
    of latitudes.
    """
    # Require odd number of latitudes
    lat = np.require(lat, dtype=np.double, requirements="C")
    assert lat.ndim == 1
    assert lat.size % 2 == 1
    # Require even number of longitudes
    lon = np.asarray(lon) # not an input to the C extension
    assert lon.ndim == 1
    assert lon.size % 2 == 0
    # Variadic function, takes multiple fields and regrids all of them
    narg = len(fields)
    assert narg > 0
    # Output shape:
    nlat_h = (lat.size - 1) // 2 + 1
    nlon_h = lon.size // 2
    # Prepare input and output fields for variadic arguments (given in
    # alternating fashion: in1, out1, in2, out2, ...)
    field_args = []
    for field in fields:
        nlev, nlat, nlon = field.shape
        assert nlat == lat.size
        assert nlon == lon.size
        field_args.append( np.require(field, dtype=np.double, requirements="C") )
        field_args.append( np.empty((nlev, nlat_h, nlon_h), dtype=np.double, order="C") )
    # Call compiled function
    status = _regrid.lib.half_resolution(
        _double_pointer(lat),
        nlev,
        nlat,
        nlon,
        narg,
        *map(_double_pointer, field_args)
    )
    if status != 0:
        raise Exception("return code {}".format(status))
    # Return half-resolution latitude and longitude arrays and all
    # half-resolution fields
    return (lat[::2], lon[::2]) + tuple(field_args[1::2])

