#include <stdio.h>
#include <stdarg.h>
#include <math.h>


#define TOL_LAT 0.001
#define DEG2RAD(angle) ((angle) / 180. * M_PI)
#define SINDEG(x) ( sin(DEG2RAD(x)) )

// Index accessors for flattened 2D fields. Can be used where the index of
// latitude is called i and the index of longitude is called j (which is the
// convention followed here).
// IDXNWW  ...  IDXNW - IDXN - IDXNE  ... IDXNEE
// IDXWW   ...  IDXW  - IDX0 - IDXE   ... IDXEE
// IDXSWW  ...  IDXSW - IDXS - IDXSE  ... IDXSEE
#define IDXNWW ((i-1) * nlon )         // Western boundary of South latitude
#define IDXNW  ((i-1) * nlon + j - 1 ) // Gridpoint Northwest
#define IDXN   ((i-1) * nlon + j )     // Gridpoint North
#define IDXNE  ((i-1) * nlon + j + 1 ) // Gridpoint Northeast
#define IDXNEE ( i    * nlon - 1 )     // Eastern boundary of North latitude
#define IDXWW  ( i * nlon )            // Western boundary of latitude
#define IDXW   ( i * nlon + j - 1 )    // Gridpoint West
#define IDX0   ( i * nlon + j )        // Gridpoint
#define IDXE   ( i * nlon + j + 1 )    // Gridpoint East
#define IDXEE  ( i * nlon + nlon - 1 ) // Eastern boundary of latitude
#define IDXSWW ((i+1) * nlon )         // Western boundary of North latitude
#define IDXSW  ((i+1) * nlon + j - 1 ) // Gridpoint Southwest
#define IDXS   ((i+1) * nlon + j )     // Gridpoint South
#define IDXSE  ((i+1) * nlon + j + 1 ) // Gridpoint Southeast
#define IDXSEE ((i+2) * nlon - 1 )     // Eastern boundary of South latitude

// Index accessor for flattened 2D fields at half resolution
#define HDX0  ((i/2)*nlon_h + (j/2))


// Sum over n elements of arr
double _sum(double * arr, size_t n) {
    double accumulator = 0.;
    for (size_t i = 0; i < n; ++i) {
        accumulator += arr[i];
    }
    return accumulator;
}


// Additional assumption: number of latitudes is odd
int half_resolution(
    double * lat,
    size_t nlev,
    size_t nlat,
    size_t nlon,
    size_t narg,
    ... // variadic
) {
    // 
    if (nlat % 2 == 0) return 1;
    if (nlon % 2 != 0) return 2;

    // Precompute area weights for gridpoint sections
    double area_north[nlat];
    double area_south[nlat];
    {
        size_t i = 0;
        area_north[i] = 0.;
        area_south[i] = SINDEG(lat[i]) - SINDEG(0.5 * (lat[i] + lat[i+1]));
        while (++i < nlat - 1) {
            area_north[i] = SINDEG(0.5 * (lat[i-1] + lat[i])) - SINDEG(lat[i]);
            area_south[i] = SINDEG(lat[i]) - SINDEG(0.5 * (lat[i] + lat[i+1]));
        }
        area_north[i] = SINDEG(0.5 * (lat[i-1] + lat[i])) - SINDEG(lat[i]);
        area_south[i] = 0.;
    }

    // Precompute normalization weights
    double normalization[nlat];
    {
        size_t i = 0;
        normalization[i] = 0.25 / (area_south[i] + area_north[i+1]);
        while (++i < nlat - 1) {
            normalization[i] = 0.25 / (area_south[i-1] + area_north[i] + area_south[i] + area_north[i+1]);
        }
        normalization[i] = 0.25 / (area_south[i-1] + area_north[i]);
    }

    // Shape of the half-resolution grid
    size_t nlat_h = 1 + (nlat - 1) / 2;
    size_t nlon_h = nlon / 2;

    // ...
    va_list args;
    va_start(args, narg);
    for (size_t k = 0; k < narg; ++k) {
        double * val_field = va_arg(args, double *);
        double * out_field = va_arg(args, double *);
        for (size_t l = 0; l < nlev; ++l) {
            double * val = val_field + (nlat   * nlon  ) * l;
            double * out = out_field + (nlat_h * nlon_h) * l;

            // Northernmost row of gridpoints:
            // ---A-----+-----B-----+---    0
            //  a | a a | b b | b b | a     → area_south[0]
            //  a | a a | b b | b b | a     → area_north[1]
            // ---+-----+-----+-----+---    1
            // .. | ... | ... | ... | ..
            {
                size_t i = 0;
                // Western (periodic) boundary
                size_t j = 0;
                out[HDX0] = normalization[i] * (
                      (val[IDXEE ] + 2.*val[IDX0] + val[IDXE ]) * area_south[i]
                    + (val[IDXSEE] + 2.*val[IDXS] + val[IDXSE]) * area_north[i+1]
                );
                // Interior
                for (size_t j = 2; j < nlon; j += 2) {
                    out[HDX0] = normalization[i] * (
                          (val[IDXW ] + 2.*val[IDX0] + val[IDXE ]) * area_south[i]
                        + (val[IDXSW] + 2.*val[IDXS] + val[IDXSE]) * area_north[i+1]
                    );
                }
                // If the northernmost row of gridpoints is the north pole,
                // enforce a single (average) value for all polar gridpoints
                if (fabs(lat[i] - 90.) < TOL_LAT) {
                    double mean_val = _sum(out+HDX0, nlon_h) / ((double) nlon_h);
                    for (size_t j = 0; j < nlon; j += 2) {
                        out[HDX0] = mean_val;
                    }
                }
            }

            // Interior:
            // .. | ... | ... | ... | ..
            // ---+-----+-----+-----+---    i-1
            //  a | a a | b b | b b | a     → area_north[i-1]
            //  a | a a | b b | b b | a     → area_south[ i ]
            // ---A-----+-----B-----+---    i
            //  a | a a | b b | b b | a     → area_north[ i ]
            //  a | a a | b b | b b | a     → area_south[i+1]
            // -- +-----+-----+-----+---    i+1
            // .. | ... | ... | ... | ..
            for (size_t i = 2; i < nlat - 1; i += 2) {
                // Western (periodic) boundary
                size_t j = 0;
                out[HDX0] = normalization[i] * (
                      (val[IDXNEE] + 2. * val[IDXN] + val[IDXNE]) *  area_south[i-1]
                    + (val[IDXEE ] + 2. * val[IDX0] + val[IDXE ]) * (area_north[i] + area_south[i])
                    + (val[IDXSEE] + 2. * val[IDXS] + val[IDXSE]) *  area_north[i+1]
                );
                // Interior
                for (size_t j = 2; j < nlon; j += 2) {
                    out[HDX0] = normalization[i] * (
                          (val[IDXNW] + 2. * val[IDXN] + val[IDXNE]) *  area_south[i-1]
                        + (val[IDXW ] + 2. * val[IDX0] + val[IDXE ]) * (area_north[i] + area_south[i])
                        + (val[IDXSW] + 2. * val[IDXS] + val[IDXSE]) *  area_north[i+1]
                    );
                }
            }

            // Southern edge
            // .. | ... | ... | ... | ..
            // ---+-----+-----+-----+---    nlat-2
            //  a | a a | b b | b b | a     → area_south[nlat-2]
            //  a | a a | b b | b b | a     → area_north[nlat-1]
            // ---A-----+-----B-----+---    nlat-1
            {
                size_t i = nlat-1;
                // Western boundary
                size_t j = 0;
                out[HDX0] = normalization[i] * (
                      (val[IDXNEE] + 2. * val[IDXN] + val[IDXNE]) * area_south[i-1]
                    + (val[IDXEE ] + 2. * val[IDX0] + val[IDXE ]) * area_north[i]
                );
                // Interior
                for (size_t j = 2; j < nlon; j += 2) {
                    out[HDX0] = normalization[i] * (
                          (val[IDXNW] + 2. * val[IDXN] + val[IDXNE]) * area_south[i-1]
                        + (val[IDXW ] + 2. * val[IDX0] + val[IDXE ]) * area_north[i]
                    );
                }
                // If the southernmost row of gridpoints is the south pole,
                // enforce a single (average) value for all polar gridpoints
                if (fabs(lat[i] + 90.) < TOL_LAT) {
                    double mean_val = _sum(out+HDX0, nlon_h) / ((double) nlon_h);
                    for (size_t j = 0; j < nlon; j += 2) {
                        out[HDX0] = mean_val;
                    }
                }
            }

        }
    }
    va_end(args);
    return 0;
}

