//
// Created by caucher on 2021/11/28.
//

#ifndef MULGIFT_RANDDATAGENERATOR_H
#define MULGIFT_RANDDATAGENERATOR_H
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <getopt.h>
#include <fstream>
#include <cstring>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define PRODUCT "TSutils - Time Series Generator\n\
Copyright (C) 2012 University of Trento\n\n"
#define STD 1   // Standard deviation
using namespace std;

class RandDataGenerator {

public:
    static void z_normalize(float *ts, int size);

    static float *generate(float *ts, int size, gsl_rng *r);

    static void generate_random_timeseries(int length, int number_of_timeseries, const char *filename);

};








#endif //MULGIFT_RANDDATAGENERATOR_H
