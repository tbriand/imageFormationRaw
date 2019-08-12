#ifndef COMPUTE_CORE
#define COMPUTE_CORE

#include <stdlib.h>
#include <math.h>
#include "fail.h"

static int bound(int x, int min, int max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

static void crop(float *out, int *cw, int *ch, float *in, int w, int h, int pd,
                 int x0, int y0, int xf, int yf)
{
    if (xf <= 0)
        xf = w + xf;
    if (yf <= 0)
        yf = h + yf;
    x0 = bound(x0, 0, w);
    xf = bound(xf, 0, w);
    y0 = bound(y0, 0, h);
    yf = bound(yf, 0, h);
    if (x0 >= xf)
        fail("bad crop x");
    if (y0 >= yf)
        fail("bad crop y");

    *cw = xf - x0;
    *ch = yf - y0;

    for (int l = 0; l < pd; l++)
        for (int j = 0; j < *ch; j++)
            for (int i = 0; i < *cw; i++)
                out[i + j*(*cw) + l*(*cw)*(*ch)] = in[i+x0 + (j+y0)*w + l*w*h];
}

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

static void affine0255(double *in, double *out, int w, int h, int pd, int type)
{
    double M = max(in, w*h);

    if( type == 0) {
        double m = min(in, w*h);
        for( int i = 0; i < w*h*pd; i++)
            out[i] = 255*(in[i]-m)/(M-m);
    }
    else {
        double factor = type*1.0/M;
        for( int i=0; i<w*h*pd; i++)
        out[i] = factor*in[i];
    }
}

#endif