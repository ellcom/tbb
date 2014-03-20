//
//  main.hpp
//  Pearson Correlation Coefficient
//  http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
//
//  Created by Elliot Adderton on 13/03/2014.
//

#ifndef A1_main_hpp
#define A1_main_hpp

#include <iostream>
#include <cmath>
#include <numeric>
#include <iomanip>

#include <tbb/tbb.h>

using namespace tbb;
using namespace std;

double mean(double*);
double standardDeviation(double *, double*);
double pearsonCorrelationCoefficient(double *, double*);


double parallelMean(double *);
double parallelStandardDeviation(double *, double *);
double parallelPearsonCorrelationCoefficient(double *, double *);

#endif
