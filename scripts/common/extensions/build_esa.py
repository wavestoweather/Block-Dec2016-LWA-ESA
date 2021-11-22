import cffi


if __name__ == "__main__":

    ffibuilder = cffi.FFI()

    code = """
        void esa_corr(double * source, double * target, size_t nmem, size_t nsrc, double * corr);
    """
    ffibuilder.cdef(code)

    ffibuilder.set_source(
        "scripts.common._esa",
        code,
        sources=[
            "scripts/common/extensions/esa.c",
        ],
        libraries=["m"],
        # Enable OpenMP support
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'],
    )
    ffibuilder.compile(verbose=True)

