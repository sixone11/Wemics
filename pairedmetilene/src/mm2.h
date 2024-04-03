#ifndef GSLMM_H
#define GSLMM_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <vector>

#define __epsilon__reserve 1e-12

/**
 *
 */
struct params{
    double * y;
    int n;
    double sy2;
};

///////// INLINE 
inline double MSIGN(double x){return(x> __epsilon__reserve ? 1.0:-1.0);}
inline double chisq_mix_P(double x, double ratio);

double sam_function(double * y, double * para, int n, double sy2);
//int gsl_sam_function(double * y, void * param, int n, double sy2, gsl_vector *f);
double gsl_sam_function(const gsl_vector * y, void * param);
double solver_sam_function(double * y, int n, double * param_init);
void solver_sam_function(double * y, int n, double * param_init, double *result, size_t control=1);
double LRF(double *x, double mu, double sigma2, size_t n);
double LRFstat(double *y, size_t n);
std::vector<double> LRFstat(std::vector<double> y);

#endif