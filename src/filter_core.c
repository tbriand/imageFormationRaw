// utility functions regarding the application of filters

/* This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) 2014-2019, Thibaud Briand <thibaud.briand@enpc.fr>
 *
 * All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "fft_core.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

static void multiply_filter(fftw_complex *outhat, fftw_complex *inhat, double *filter, int w, int h, int pd)
{    
    double tmp;
    for (int i = 0; i < w*h; i++) {
            tmp = filter[i];
            for (int l = 0; l < pd; l++)
                outhat[i + l*w*h] = inhat[i + l*w*h]*tmp;
    }
}

static void compute_asymptotic_nc_filter(double *filter, int order, double sigma, int w, int h)
{
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    double fw = 2*M_PI/((double) w);
    double fh = 2*M_PI/((double) h);
    double sigma2 = sigma*sigma*0.5;
    double x, y;
    
    if ( order <= 1 ) {
        for (int j = 0; j < h; j++) {
            y = (j<ind_h) ? j : j-h;
            y *= fh;
            for (int i = 0; i < w; i++) {
                x = (i<ind_w) ? i: i-w;
                x *= fw;
                x = x*x+y*y;
                filter[i+j*w] = exp(-sigma2*x);
            }
        }
    }
    else { //order = 2
        for (int j = 0; j < h; j++) {
            y = (j<ind_h) ? j : j-h;
            y *= fh;
            for (int i = 0; i < w; i++) {
                x = (i<ind_w) ? i: i-w;
                x *= fw;
                x = x*x+y*y;
                filter[i+j*w] = exp(-sigma2*x)*(1+sigma2*x);
            }
        }
    }
}

static void apply_asymptotic_nc_filter(double *in, double *out, int w, int h, int pd, int order, double sigma, int inverse)
{
    fftw_complex *inhat = fftw_malloc(w*h*pd*sizeof*inhat);
    fftw_complex *outhat = fftw_malloc(w*h*pd*sizeof*outhat);
    double *filter = malloc(w*h*sizeof*filter);
    
    // compute input DFT
    do_fft_real(inhat, in, w, h, pd);

    // compute filter
    compute_asymptotic_nc_filter(filter, order, sigma, w, h);
    
    // inverse if needed
    if( inverse == 1 ) {
        for (int i=0; i<w*h; i++)
            filter[i] = 1.0/filter[i];
    }
        
    // multiply filter
    multiply_filter(outhat, inhat, filter, w, h, pd);

    // DFT inverse
    do_ifft_real(out, outhat, w, h, pd);
    
    // free memory
    fftw_free(inhat);
    fftw_free(outhat);
}

static void symmetrization(double *in, double *out, int w, int h, int pd)
{
    for(int l=0; l<pd; l++) {
        /* Top left corner */
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++)
                out[i+j*2*w+l*4*w*h] = in[i+j*w+l*w*h];
    
        /* Top right corner */
        for(int j=0; j<h; j++)
            for(int i=w; i<2*w; i++)
                out[i+j*2*w+l*4*w*h] = in[2*w-1-i + j*w + l*w*h];
        
        /* Bottom left corner */
        for(int i=0; i<w; i++)
            for(int j=h; j<2*h;j++)
            out[i+j*2*w+l*4*w*h] = in[i + (2*h-1-j)*w + l*w*h];

        /* Bottom right corner */
        for(int i=w; i<2*w; i++)
            for(int j=h; j<2*h;j++)
                out[i+j*2*w+l*4*w*h] = in[2*w-1-i + (2*h-1-j)*w + l*w*h];
    }
}

void apply_asymptotic_nc_filter_global(double *in, double *out, int w, int h, int pd, int order, double sigma, int inverse)
{
    // symmetrization of the input
    double *in_sym = malloc(4*w*h*pd*sizeof*in_sym);  
    symmetrization(in, in_sym, w, h, pd);
    
    // apply filter to the symmetrized image
    apply_asymptotic_nc_filter(in_sym, in_sym, 2*w, 2*h, pd, order, sigma, inverse);
    
    // keep part of the image
    for(int l = 0; l < pd; l++)
        for(int i = 0; i < w; i++)
            for(int j = 0; j < h; j++)
                out[i + j*w + l*w*h] = in_sym[i + j*2*w + l*4*w*h];

    // free memory
    free(in_sym);
}

static void bilinear_demosaicing(double *in, double *out, int w, int h) {
    //filters for bilinear interp
    int RB[3][3] = {
        {1,2,1},
        {2,4,2},
        {1,2,1}
    };
    int G[3][3] = {
        {0,1,0},
        {1,4,1},
        {0,1,0}
    };

    //CFA image with 0 when unknown
    double *CFA = malloc(3*w*h*sizeof*CFA);
    for(int i = 0; i < 3*w*h; i++)
        CFA[i] = 0;
    
    int l;
    for(int j = 0; j < h; j++)
        for(int i = 0; i < w; i++) {
            l = i%2+j%2;
            CFA[i + j*w + l*w*h] = in[j*w+i];
        }

    //bilinear interp
    double tmpR, tmpG, tmpB;
    for(int j = 0; j < h; j++)
        for(int i = 0; i < w; i++) {
            tmpR = 0;
            tmpG = 0;
            tmpB = 0;
            for (int x=-1; x<=1 ; x++) {
                int indx=j+x;
                for (int y=-1; y<=1 ; y++) {
                    int indy=i+y;
                    if( indx>=0 && indx<h && indy>=0 && indy<w ) {
                        tmpR += CFA[indy+indx*w]*RB[y+1][x+1];
                        tmpG += CFA[indy+indx*w+w*h]*G[y+1][x+1];
                        tmpB += CFA[indy+indx*w+2*w*h]*RB[y+1][x+1];
                    }
                }
            }
            out[i+j*w] = tmpR*0.25;
            out[i+j*w+w*h] = tmpG*0.25;
            out[i+j*w+2*w*h] = tmpB*0.25;
        }
    
    //free memory
    free(CFA);
}

void apply_bilinear_luminance(double *in, double *out, int w, int h) {
    // memory allocation
    double *rgb = malloc(3*w*h*sizeof*rgb);
    
    // bilinear demosaicing
    bilinear_demosaicing(in, rgb, w, h);
    
    // average of the channels
    for(int i=0; i<w*h; i++)
        out[i] = 0.3333333333333*(rgb[i] + rgb[i+w*h] + rgb[i+2*w*h]);
    
    // free memory
    free(rgb);
}

void apply_hvs_lowpass(double *double_in, double *double_out, int w, int h)
{
    int filter[11][11] = {
        {0,0,0,0,1,0,1,0,0,0,0},
        {0,0,0,-1,0,-2,0,-1,0,0,0},
        {0,0,1,1,2,1,2,1,1,0,0},
        {0,-1,1,-5,3,-9,3,-5,1,-1,0},
        {1,0,2,3,1,7,1,3,2,0,1},
        {0,-2,1,-9,7,104,7,-9,1,-2,0},
        {1,0,2,3,1,7,1,3,2,0,1},
        {0,-1,1,-5,3,-9,3,-5,1,-1,0},
        {0,0,1,1,2,1,2,1,1,0,0},
        {0,0,0,-1,0,-2,0,-1,0,0,0},
        {0,0,0,0,1,0,1,0,0,0,0}
    };

    double ifactor = 1.0/128;
    double tmp;
    for(int j = 0; j < h; j++)
        for(int i = 0; i < w; i++) {
            tmp=0.0;
            for (int x=-5; x<=5 ; x++) {
                int indx=j+x;
                for (int y=-5; y<=5 ; y++) {
                    int indy=i+y;
                    if( indx>=0 && indx<h && indy>=0 && indy<w )
                        tmp += double_in[indy+indx*w]*filter[x+5][y+5];
                }
            }
            double_out[i+j*w] = tmp*ifactor;
        }
}

static void lowpass(fftw_complex *inhat, fftw_complex *outhat, int w, int h, int pd)
{
    // fftshift
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    int i3, j3;

    for(int l=0; l<pd; l++)
        for(int j=0; j<h; j++) {
            j3 = (j<ind_h) ? j : j - h;
            for(int i=0; i<w; i++) {
                i3 = (i<ind_w) ? i: i - w;
                outhat[i+j*w+l*w*h]= ( abs(i3) < ind_w/2 && abs(j3) < ind_h/2 ) ? inhat[i+j*w+l*w*h] : 0;
            }
        }
}

void apply_perfect_lowpass(double *in, double *out, int w, int h, int pd)
{
    // memory allocation
    fftw_complex *inhat = malloc(w*h*pd*sizeof*inhat);
    fftw_complex *outhat = malloc(w*h*pd*sizeof*outhat);

    // compute DFT of the input
    do_fft_real(inhat, in, w, h, pd);

    // set to zero the values */
    lowpass(inhat, outhat, w, h, pd);

    // compute iDFT of the input
    do_ifft_real(out, outhat, w, h, pd);

    // free memory
    fftw_free(inhat);
    fftw_free(outhat);
}

#define sinc(x) x ? sin(x)/(x) : 1
static void integration_kernel(fftw_complex *inhat, fftw_complex *outhat, int w, int h, int pd, double zoom)
{
    // fftshift
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    int i3, j3;
    double factori, factorj;
    
    for(int l=0; l<pd; l++)
        for(int j=0; j<h; j++) {
            j3 = (j<ind_h) ? j : j - h;
            factorj = sinc(zoom*M_PI/h*j3);
            for(int i=0; i<w; i++) {
                i3 = (i<ind_w) ? i: i - w;
                factori = sinc(zoom*M_PI/w*i3);
                outhat[i+j*w+l*w*h] = factori*factorj*inhat[i+j*w+l*w*h];
            }
        }
}

void apply_integration_kernel(double *in, double *out, int w, int h, int pd, double zoom)
{
    // memory allocation
    fftw_complex *inhat = malloc(w*h*pd*sizeof*inhat);
    fftw_complex *outhat = malloc(w*h*pd*sizeof*outhat);

    // compute DFT of the input
    do_fft_real(inhat, in, w, h, pd);

    // apply the integration kernel */
    integration_kernel(inhat, outhat, w, h, pd, zoom);

    // compute iDFT of the input
    do_ifft_real(out, outhat, w, h, pd);

    // free memory
    fftw_free(inhat);
    fftw_free(outhat);
}