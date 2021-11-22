import numpy as np
import cffi
from . import _esa


# CFFI preparation
_ffi = cffi.FFI()

def _double_pointer(arr):
    return _ffi.cast("double *", _ffi.from_buffer(arr))


# Bindings to compiled correlation-based ESA
def esa_corr(src, tgt):
    nmem, *shp = src.shape
    src = src.reshape((nmem, -1))
    nfld = np.product(shp)
    assert nmem == tgt.size
    source = np.require(src, dtype=np.double, requirements="C")
    target = np.require(tgt, dtype=np.double, requirements="C")
    output = np.empty(nfld, dtype=np.double, order="C")
    _esa.lib.esa_corr(
        _double_pointer(source),
        _double_pointer(target),
        nmem,
        nfld,
        _double_pointer(output)
    )
    return output.reshape(shp)

def esa_corr_with_confidence(src, tgt, nsamples=1000, sample_size=None):
    if sample_size is None:
        sample_size = tgt.shape[0]
    # Need at least two samples to obtain a confidence estimate
    if nsamples < 2:
        return esa_corr(src, tgt), None
    # Bootstrap sampling
    rng = np.random.default_rng()
    samples = []
    for _ in range(nsamples):
        i = rng.choice(np.arange(tgt.shape[0]), size=sample_size, replace=True)
        samples.append(esa_corr(src[i,...], tgt[i,...]))
    samples = np.asarray(samples)
    # Confidence interval given by one standard deviation of the correlation
    # coeffienct distribution
    return esa_corr(src, tgt), np.std(samples, axis=0)

