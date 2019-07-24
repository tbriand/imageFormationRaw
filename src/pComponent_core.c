#ifndef PCOMPONENT_CORE
#define PCOMPONENT_CORE

#include "fft_core.c"

static void jumps(double *out, double *in, int w, int h, int pd)
{
    double (*x)[w][pd] = (void*)in;
    double (*y)[w][pd] = (void*)out;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
            for (int l = 0; l < pd; l++)
                y[j][i][l]=0;
    for (int j = 0; j < h; j++)
        for (int l = 0; l < pd; l++)
        {
            y[j][0][l] = x[j][0][l] - x[j][w-1][l];
            y[j][w-1][l] = -y[j][0][l];
        }
    for (int i = 0; i < w; i++)
        for (int l = 0; l < pd; l++)
        {
            y[0][i][l] += x[0][i][l] - x[h-1][i][l];
            y[h-1][i][l] -= x[0][i][l] - x[h-1][i][l];
        }
}

static void compute_sComponenthat(double *in, double *shat, int w, int h, int pd)
{
    double *v = xmalloc(w*h*pd*sizeof*v);
    double *vhat = xmalloc(w*h*2*pd*sizeof*vhat);
    jumps(v, in, w, h, pd);
    double (*x)[w][2*pd] = (void*)vhat;
    double (*y)[w][2*pd] = (void*)shat;
    fft_direct(vhat, v, w, h, pd);
    
    double tmp;
    double factorh = 2*M_PI/h;
    double factorw = 2*M_PI/w;
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            tmp = 1.0/(4-2*cos(j*factorh)-2*cos(i*factorw));
            for (int l = 0; l < 2*pd; l++)
                y[j][i][l] = x[j][i][l]*tmp;
        }
    normalize_double_array_inplace(shat, w*h*2*pd);

    free(v);
    free(vhat);
}

static void compute_pComponenthat(double *in, double *shat, double *phat, int w, int h, int pd)
{
    double *inhat = xmalloc(w*h*2*pd*sizeof*inhat);
    fft_direct(inhat, in, w, h, pd);
    for (int i=0; i<2*pd*w*h; i++)
    phat[i] = inhat[i]-shat[i];
    free(inhat);
}

void compute_component(double *in, double *pComponent_zoom, double *sComponent, int zoom, int w, int h, int pd)
{
    double *shat = xmalloc(w*h*2*pd*sizeof*shat);
    double *phat = xmalloc(w*h*2*pd*sizeof*phat);
    double *phat_zero = xmalloc(zoom*zoom*w*h*2*pd*sizeof*phat_zero);
    
    // compute smooth component
    compute_sComponenthat(in, shat, w, h, pd);
    fft_inverse(sComponent, shat, w, h, 2*pd);

    // compute the zoomed periodic component
    // 1) fft image - fft sComponent
    compute_pComponenthat(in, shat, phat, w, h, pd);
    // 2) zero-padding
    zoom_in(phat, phat_zero, zoom, w, h, pd, 1);
    // 3) fft inverse
    fft_inverse(pComponent_zoom, phat_zero, zoom*w, zoom*h, 2*pd);

    free(shat);
    free(phat);
    free(phat_zero);
}

#endif

