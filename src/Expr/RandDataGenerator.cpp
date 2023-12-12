//
// Created by caucher on 2021/11/28.
//

#include "../../include/Expr/RandDataGenerator.h"
void RandDataGenerator::z_normalize(float *ts, int size) {
    int i;
    float mean = 0;//gsl_stats_mean(ts, 1, size);
    float std = 0;//gsl_stats_sd(ts, 1, size);
    for (i = 0; i < size; i++) {
        mean += ts[i];
    }
    mean /= size;

    for (i = 0; i < size; i++) {
        std += (ts[i] - mean) * (ts[i] - mean);
    }
    std /= size;
    std = sqrt(std);
    for (i = 0; i < size; i++) {
        ts[i] = (ts[i] - mean) / std;
    }
}

float *RandDataGenerator::generate(float *ts, int size, gsl_rng *r) {
    int i;
    float x = 0, dx;

    for (i = 0; i < size; i++) {
        dx = gsl_ran_gaussian(r, STD); // mean=0, std=STD
        x += dx;
        ts[i] = x;
    }

    z_normalize(ts, size);
    return ts;
}


/**
    Generates a set of random time series.
**/
void RandDataGenerator::generate_random_timeseries(int length, int number_of_timeseries, const char *filename) {
    // Initialize random number generation
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    FILE *data_file;
    data_file = fopen(filename, "wb");

    auto *ts = new float[length];
    int i;
    for (i = 1; i <= number_of_timeseries; i++) {
        generate(ts, length, r);
        fwrite(ts, sizeof(float), length, data_file);
    }

    if (i % (1000) == 0) {
        fprintf(stderr, "\r\x1b[m>> Generating: \x1b[36m%2.2lf%%\x1b[0m",
                (float) ((float) i / (float) number_of_timeseries) * 100);
    }
    // Finalize random number generator
    fclose (data_file);
    gsl_rng_free (r);
}