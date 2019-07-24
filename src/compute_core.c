#ifndef COMPUTE_CORE
#define COMPUTE_CORE

#include <stdlib.h>
#include <math.h>

static double mean(double *in, int N) {
    
    double val = 0.0;
    for (int i=0; i<N; i++)
        val += in[i];
    val /= (double) N;
    
    return val;
}

static double rmse(double *in, int N) {
    
    double val = 0.0;
    for (int i=0; i<N; i++)
        val += in[i]*in[i];
    val /= (double) N;
    val = sqrt(val);
    
    return val;
}

static double std(double *in, double m, int N) {
    
    double val = 0.0;
    double tmp;
    for (int i=0; i<N; i++) {
        tmp = in[i] - m;
        val += tmp*tmp;
    }
    val /= (double) N;
    val = sqrt(val);
    
    return val;
}

static double max(double *in, int N) {
    
    double val = in[0];
    for (int i=1; i<N; i++)
        val = (in[i] > val) ? in[i] : val; 
    
    return val;
}

static double min(double *in, int N) {
    
    double val = in[0];
    for (int i=1; i<N; i++)
        val = (in[i] < val) ? in[i] : val; 
    
    return val;
}

static double correlation(double *im1, double *im2, int N) {
    //mean computations
    double m1 = mean(im1, N);
    double m2 = mean(im2, N);
    
    //std computations
    double std1 = std(im1, m1, N);
    double std2 = std(im2, m2, N);
    
    //correlation computation
    double val = 0.0;
    for(int i=0; i<N; i++)
        val += (im1[i] - m1)*(im2[i] - m2);
    val /= (N*std1*std2);
    
    return val;
}

#endif