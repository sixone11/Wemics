#include <stdio.h>
#include <math.h>
#include "mm2.h"  //add by xp
#include <iostream>  //add by xp
#include <gsl/gsl_cdf.h>  //add by xp


#define RETRY 100
#define MIN_ITER 50
#define MAX_ITER 1000

double sam_function(double * y, double * para, int n, double sy2){
    double s = 0;
    double mu_pi = 3.1415926535897932;
    size_t i = 0;
    // n = len(y)
    for (;i<n; i++){
        s = s + log(cosh(para[0] * y[i] / para[1]));
    }
    s = -((double) n)* log(2.0/mu_pi/para[1])/2. +
        (double) n * para[0] * para[0] / 2.0/para[1] + sy2/2.0/para[1] -
        s;
    return s;
}

// In GSL, param is the constant variable one uses in the likelihood
// const gsl_vector * y, meanwhile, is your parameters to be updated ...
double gsl_sam_function(const gsl_vector * y, void * param){
    double para[2] = {0.,0.};
    para[0] = gsl_vector_get(y,0); //mean
    para[1] = gsl_vector_get(y,1); //sd
    return(
        sam_function(
                ( (struct params *) param)->y,
                para,
                ( (struct params *) param)->n,
                ( (struct params *) param)->sy2)
        );
}

void solver_sam_function(double * y, int n, double * param_init, double *result, size_t control){
    gsl_set_error_handler_off();
    double sy2=0;
    for(size_t ii =0;ii<n;ii++) sy2+=y[ii]*y[ii];
    size_t iter=0;
    double size=0;
    int status;
    size_t status_loop = 0;
    struct params p = {y, n, sy2};
    gsl_multimin_function f = {&gsl_sam_function, 2, &p};
    gsl_vector *ss, *x;
    double x_init[2] = {param_init[0], param_init[1]};
    x = gsl_vector_alloc(2);
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2rand;//gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    // Initialization
    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);
    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 1.0);

    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &f, x, ss);

        try{
            do{ status_loop=0;
                iter++;
                status = gsl_multimin_fminimizer_iterate(s);

                if (status && iter <= MIN_ITER){
                    break;
                }

                size = gsl_multimin_fminimizer_size (s);
                status = gsl_multimin_test_size (size, 1e-2);
                if (control==1){
                    if (status == GSL_SUCCESS) printf ("converged to minimum at\n");
                    printf ("%10d %10.6e %10.6e f() = %10.6f size = %.6f\n",
                            (int)iter,
                            gsl_vector_get (s->x, 0),
                            gsl_vector_get (s->x, 1),
                            s->fval, size);
                }
                if (iter <= MIN_ITER) status_loop = 1;
                if (iter < MAX_ITER && status == GSL_CONTINUE) status_loop = 1;
                } while ( status_loop == 1);
        }
        catch (...){
            printf ("Infinite");
        }

  result[0] = gsl_vector_get (s->x, 0);
  result[1] = gsl_vector_get (s->x, 1);
  result[2] = (double) status;

    gsl_vector_free (x);
    gsl_vector_free (ss);
    gsl_multimin_fminimizer_free (s);
}

double LRF(double *x, double mu, double sigma2, size_t n){
    double pi=3.1415926535897932;
    double re0=0.,re1=0.,re=0;
    size_t i = 0;
    if (sigma2<1e-14) {
        sigma2=1;
    }
    for(;i<n;i++){
        re0 = re0 + x[i] * x[i];
        re1 = re1 + log(cosh(mu*x[i]/sigma2));
    }
    re = (double)n/2.*log(2/pi) -  (double)n/2.*log(sigma2) - (re0 + n*mu*mu)/sigma2/2 + re1;
return re;
}

/*
*  -2 Log-likelihood
*/
double LRFstat(double* y, size_t n) {
    double mean = 0.;
    double sd = 0.;
    double mle[3] = { 0., 0., 0. };
    double lrt;
    for (size_t i = 0; i < n; i++) mean += y[i];
    mean = mean / n;
    for (size_t i = 0; i < n; i++) sd += (y[i] - mean) * (y[i] - mean);
    sd = sqrt(sd / (n + 1.));
    double param_init[2] = { mean, sd };
    try {
        solver_sam_function(y, n, param_init, mle, 0);
    }
    catch (...) {
    }
    if (mle[2] > -0.5) {
        double meany2 = 0;
        for (size_t i = 0; i < n; i++) {
            meany2 += y[i] * y[i];
        }
        meany2 = meany2 / n;
        lrt = 2 * (LRF(y, mle[0], mle[1], n) - LRF(y, 0., meany2, n));
        lrt = lrt < 0 ? 0 : lrt;
    }
    else {
        lrt = -1.;
    }
   // return lrt;
    double chisq_Q = 1 - chisq_mix_P(lrt, 0.5);
    return chisq_Q;
}

/*
*  -2 Log-likelihood remask
*/
std::vector<double> LRFstat(std::vector<double> input_array){
    std::vector<double> result(3,0.);
    double * y = &input_array[0];
    int n = input_array.size();
        double mean = 0.;
    double sd = 0.;
    double mle[3] = {0., 0., 0.};
    double lrt;
    for(size_t i = 0; i<n; i++) mean += y[i];
    mean = mean/n;
    for(size_t i = 0; i<n; i++) sd += (y[i]-mean)*(y[i]-mean);
    sd = sqrt(sd/(n+1.));
    double param_init[2] = {mean, sd};
        try{
            solver_sam_function(y, n, param_init, mle, 0);
        }catch(...){
        }
    if(mle[2] > -0.5){
        double meany2 = 0;
        for (size_t i = 0; i<n;i++){
            meany2 += y[i]*y[i];
        }
        meany2 = meany2/n;
        lrt = 2 *(LRF(y, mle[0], mle[1], n)-LRF(y, 0., meany2, n));
        lrt= lrt<0?0:lrt;
    }else{
        lrt = -1.;
    }
    result[0] = mle[0];
    result[1] = mle[1];
    result[2] = 1 - chisq_mix_P(lrt, 0.5);
    return result;
}

/*
f(x) MUST be monotonic
*/
double f_uniroot(double f(double,std::vector<double>), double start, double end, std::vector<double> paramas, bool debug){
    // size_t MAX_ITER = 2000;
    double split_lambda = 0.618;
    double split_lambda2 = 0.382;
    double tol = 1e-12;
    double f_l, f_r, f_l0;
    size_t iter =0;
    double tmp = end; 
    if (end < start) {end = start; start = tmp;}
    f_l = f(start, paramas);
    f_r = f(end, paramas);
    do{        
        if ( (MSIGN(f_l) ) * (MSIGN(f_r)) > 0.5) break;
        tmp = split_lambda*start + split_lambda2*end;
        f_l0 = f(tmp, paramas);
        if (fabs(f_l0) < tol) {iter = 10000; break;}
        if ( ( MSIGN(f_l0))  * (MSIGN(f_r)) > 0.5) {
            end = tmp;
            f_r = f_l0;
            if (debug) printf("+\t");
        }else{
            start = tmp;
            f_l = f_l0;
            if (debug) printf("-\t");
        }
        if (debug) printf("Iter %3d f(%4f,%4f)=(%4f, %4f).\n", (int) iter, start, end, f_l, f_r);
        iter++;
    }while (iter < MAX_ITER);
    if (iter < 2) printf("Same Sign, choose another initial start point.\n");
    return(end);
}

inline double chisq_mix_P(double x, double ratio){
    x = (x < 0 ? 0:x);
    ratio = (ratio > 1 ? 1:ratio);
    ratio = (ratio < 0 ? 0:ratio);
    return( x < 0 ? 0: (ratio * gsl_cdf_chisq_P(x, 0) + (1 - ratio) * gsl_cdf_chisq_P(x, 1)) );
}