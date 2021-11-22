import cffi


if __name__ == "__main__":

    ffibuilder = cffi.FFI()

    code = """
        int half_resolution(double * latitude, size_t nlev, size_t nlat, size_t nlon, size_t narg, ...);
    """
    ffibuilder.cdef(code)

    ffibuilder.set_source(
        "scripts.common._regrid",
        code,
        sources=[
            "scripts/common/extensions/regrid.c",
        ],
        libraries=["m"]
    )
    ffibuilder.compile(verbose=True)

