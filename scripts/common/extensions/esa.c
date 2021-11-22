#include <stdio.h>
#include <math.h>


void esa_corr(double * source, double * target, size_t nmem, size_t nsrc, double * corr) {

    // Mean value of target
    double target_avg = 0.;
    for (size_t j = 0; j < nmem; ++j) {
        target_avg += target[j];
    }
    target_avg /= (double) nmem;

    // Standard deviation of target
    double target_term[nmem];
    double target_std = 0.;
    for (size_t j = 0; j < nmem; ++j) {
        target_term[j] = target[j] - target_avg;
        target_std += target_term[j] * target_term[j];
    }
    // The normalizations for std and cov cancel in the correlation
    // coefficient, omit them
    target_std = sqrt(target_std);

    #pragma omp parallel for
    for (size_t i = 0; i < nsrc; ++i) {
        // Mean value of source
        double source_avg = 0.;
        for (size_t j = 0; j < nmem; ++j) {
            source_avg += source[j * nsrc + i];
        }
        source_avg /= (double) nmem;

        // Covariance and variance of source
        double covariance = 0.;
        double source_var = 0.;
        for (size_t j = 0; j < nmem; ++j) {
            double source_term = source[j * nsrc + i] - source_avg;
            covariance += source_term * target_term[j];
            source_var += source_term * source_term;
        }

        // Correlation coefficient
        corr[i] = covariance / target_std / sqrt(source_var);
    }
}
