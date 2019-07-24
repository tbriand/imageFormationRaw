// Utility functions regarding the FFT

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

#ifndef FFT_CORE
#define FFT_CORE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <fftw3.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846 /* pi */
#endif

#define FORI(n) for(int i=0;i<(n);i++)
#define FORJ(n) for(int j=0;j<(n);j++)
#define FORL(n) for(int l=0;l<(n);l++)

#ifndef USE_WISDOM
void evoke_wisdom(void) {}
void bequeath_wisdom(void) {}
// #else//USE_WISDOM
// #include "fftwisdom.c"
#endif//USE_WISDOM

#include "xmalloc.c"

// wrapper around FFTW3 that computes the complex-valued Fourier transform
// of a real-valued image
static void fft_2ddouble(fftw_complex *fx, double *x, int w, int h)
{
    fftw_complex *a = fftw_malloc(w*h*sizeof*a);

    //fprintf(stderr, "planning...\n");
    evoke_wisdom();
    fftw_plan p = fftw_plan_dft_2d(h, w, a, fx,
                                     FFTW_FORWARD, FFTW_ESTIMATE);
    bequeath_wisdom();
    //fprintf(stderr, "...planned!\n");

    FORI(w*h) a[i] = x[i]; // complex assignment!
    fftw_execute(p);

    fftw_destroy_plan(p);
    fftw_free(a);
    fftw_cleanup();
}


// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2ddouble(double *ifx,  fftw_complex *fx, int w, int h)
{
    fftw_complex *a = fftw_malloc(w*h*sizeof*a);
    fftw_complex *b = fftw_malloc(w*h*sizeof*b);

    //fprintf(stderr, "planning...\n");
    evoke_wisdom();
    fftw_plan p = fftw_plan_dft_2d(h, w, a, b,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    bequeath_wisdom();
    //fprintf(stderr, "...planned!\n");

    FORI(w*h) a[i] = fx[i];
    fftw_execute(p);
    double scale = 1.0/((double) w*h);
    FORI(w*h) {
        fftw_complex z = b[i] * scale;
        ifx[i] = creal(z);
        //assert(cimagf(z) < 0.001);
    }
    fftw_destroy_plan(p);
    fftw_free(a);
    fftw_free(b);
    fftw_cleanup();
}

// Wrapper around FFTW3 that computes the real-valued inverse Fourier transform
// of a complex-valued frequantial image.
// The input data must be hermitic.
static void ifft_2ddouble_unscaled(double *ifx,  fftw_complex *fx, int w, int h)
{
    fftw_complex *a = fftw_malloc(w*h*sizeof*a);
    fftw_complex *b = fftw_malloc(w*h*sizeof*b);

    //fprintf(stderr, "planning...\n");
    evoke_wisdom();
    fftw_plan p = fftw_plan_dft_2d(h, w, a, b,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    bequeath_wisdom();
    //fprintf(stderr, "...planned!\n");

    FORI(w*h) a[i] = fx[i];
    fftw_execute(p);
    FORI(w*h) {
        ifx[i] = creal(b[i]);
        //assert(cimagf(z) < 0.001);
    }
    fftw_destroy_plan(p);
    fftw_free(a);
    fftw_free(b);
    fftw_cleanup();
}

static void ifft_unscaled_float(float *ifx, fftw_complex *fx, int w, int h, int pd)
{
    fftw_complex *a = fftw_malloc(w*h*sizeof*a);
    fftw_complex *b = fftw_malloc(w*h*sizeof*b);

    //fprintf(stderr, "planning...\n");
    evoke_wisdom();
    fftw_plan p = fftw_plan_dft_2d(h, w, a, b,
                                     FFTW_BACKWARD, FFTW_ESTIMATE);
    bequeath_wisdom();
    //fprintf(stderr, "...planned!\n");

    FORL(pd) {
        FORI(w*h) a[i] = fx[i*pd+l];
        fftw_execute(p);
        FORI(w*h) {
            ifx[i*pd+l] = creal(b[i]);
            //assert(cimagf(z) < 0.001);
        }
    }
    
    fftw_destroy_plan(p);
    fftw_free(a);
    fftw_free(b);
    fftw_cleanup();
}

// Wrapper around FFTW3 that computes the iDCT
// of a real-valued image.
static void ifct_2ddouble(double *ifx,  double *fx, int w, int h)
{
    evoke_wisdom();
    fftw_plan p = fftw_plan_r2r_2d(h, w, fx, ifx, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    bequeath_wisdom();
    
    fftw_execute(p);
    fftw_destroy_plan(p);
    fftw_cleanup();
}



// if it finds any strange number, sets it to zero
static void normalize_double_array_inplace(double *x, int n)
{
    for (int i = 0; i < n; i++)
        if (!isnormal(x[i]))
            x[i] = 0;
}

static void fft_direct(double *y, double *x, int w, int h, int pd)
{
    double *c = xmalloc(w*h*sizeof*c);
    fftw_complex *gc = xmalloc(w*h*sizeof*gc);
    FORL(pd) {
        FORI(w*h)
            c[i] = x[i*pd + l];
        fft_2ddouble(gc, c, w, h);
        FORI(w*h) {
            y[2*(i*pd + l)+0] = creal(gc[i]);
            y[2*(i*pd + l)+1] = cimag(gc[i]);
        }
    }
    free(c);
    free(gc);
}

static void fft_inverse(double *y, double *x, int w, int h, int pd)
{
    int pdh = pd/2;
    assert(pd == 2*pdh);
    fftw_complex *c = xmalloc(w*h*sizeof*c);
    double *gc = xmalloc(w*h*sizeof*gc);
    FORL(pdh) {
        FORI(w*h)
            c[i] = x[i*pd + 2*l] + I * x[i*pd+2*l+1];
        ifft_2ddouble(gc, c, w, h);
        FORI(w*h)
            y[i*pdh + l] = gc[i];
    }
    free(c);
    free(gc);
}

static void fftshift( fftw_complex *fhat, fftw_complex *fshift, int w, int h, int pd) {
    int w_test = w/2;
    int h_test = h/2;
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    int i2, j2;
  
    for(int i=0; i<w; i++) {
        i2 = (i<ind_w) ? i + w_test : i - ind_w;
        for(int j=0; j<h; j++) {
            j2 = (j<ind_h) ? j + h_test : j - ind_h;
            for(int l=0; l<pd; l++)
                fshift[(j2*w+i2)*pd + l] = fhat[(j*w+i)*pd + l];
        }
    }
}

static void fftshift_double( double *f, double *fshift, int w, int h, int pd) {
    int w_test = w/2;
    int h_test = h/2;
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    int i2, j2;
  
    for(int i=0; i<w; i++) {
        i2 = (i<ind_w) ? i + w_test : i - ind_w;
        for(int j=0; j<h; j++) {
            j2 = (j<ind_h) ? j + h_test : j - ind_h;
            for(int l=0; l<pd; l++)
                fshift[(j2*w+i2)*pd + l] = f[(j*w+i)*pd + l];
        }
    }
}

static void ifftshift( fftw_complex *fhat, fftw_complex *fshift, int w, int h, int pd) {
    int indx, indy;
    int w_test = w/2;
    int h_test = h/2;
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    for (int i=0; i<w; i++) {
        for( int j=0; j<h; j++) {
            indx = (i < w_test ) ? i+ind_w : i-ind_w;
            indy = (j < h_test ) ? j+ind_h : j-ind_h;
            for(int l=0; l<pd; l++)
                fshift[(j*w+i)*pd+l] = fhat[(indy*w+indx)*pd+l];
        }
    }
}

static void symmetrization(double *in, double *out, int w, int h, int pd)
{
  double (*in2)[w][pd] = (void*) in;
  double (*out2)[2*w][pd] = (void*) out;

  FORL(pd){
    /* Top left corner */
    for(int i=0;i<w;i++){
      for(int j=0; j<h;j++){
	out2[j][i][l] = in2[j][i][l];
      }
    }
    /* Top right corner */
    for(int i=w;i<2*w;i++){
      for(int j=0; j<h;j++){
	out2[j][i][l] = in2[j][2*w-1-i][l];
      }
    }
    /* Bottom left corner */
    for(int i=0;i<w;i++){
      for(int j=h; j<2*h;j++){
	out2[j][i][l] = in2[2*h-1-j][i][l];
      }
    }
    /* Bottom right corner */
    for(int i=w;i<2*w;i++){
      for(int j=h; j<2*h;j++){
	out2[j][i][l] = in2[2*h-1-j][2*w-1-i][l];
      }
    }
  }
}

static void lowpass_dft(double *inhat, double *outhat, int w, int h, int pd)
{
  double (*outhat2)[w][2*pd] = (void*) outhat;
  double (*inhat2)[w][2*pd]  = (void*) inhat;

  /* fftshift */
  int ind_w = (w+1)/2;
  int ind_h = (h+1)/2;
  int i3, j3;

  FORI(w){
    i3 = (i<ind_w) ? i: i - w;
    FORJ(h){
      j3 = (j<ind_h) ? j : j - h;
      FORL(2*pd){
        outhat2[j][i][l]= ( abs(i3) < ind_w/2 && abs(j3) < ind_h/2 ) ? inhat2[j][i][l] : 0;
          }
      }
  }
}

static void apply_perfect_lowpass(double *in, double *out, int w, int h, int pd)
{
  double *inhat = xmalloc(w*h*2*pd*sizeof*inhat);
  double *outhat = xmalloc(w*h*2*pd*sizeof*outhat);

  /* compute the FFT of the input */
  fft_direct(inhat, in, w, h, pd);

  /* set to zero the values */
  lowpass_dft(inhat, outhat, w , h ,pd);

  /* compute the iFFT */
  fft_inverse(out, outhat, w, h, 2*pd);

  /* free memory */
  free(inhat);
  free(outhat);
}

static void zoom_in(double *inhat, double *outhat, int zoom, int w, int h, int pd, int option)
{
  if( zoom > 1 ) {
  double norm = zoom*zoom;
  double (*outhat2)[w*zoom][2*pd] = (void*) outhat;
  double (*inhat2)[w][2*pd] = (void*) inhat;

    int w_test = w/2;
    int h_test = h/2;
    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;

    /* Fill the output array with zeros */
     for (int i=0; i<2*pd*zoom*zoom*w*h; i++)
        outhat[i]=0.0;

    int jj, ii; /* new coordinates for original pixel (j,i) */
    /* fill the corners with the values */
    for(int j=0; j< h; j++){
        jj = (j<ind_h) ? j : j + (zoom-1)*h;
        for(int i=0; i<w; i++){
            ii = (i<ind_w) ? i : i + (zoom-1)*w;
            for (int l = 0; l < 2*pd; l++)
                outhat2[jj][ii][l]  = inhat2[j][i][l]*norm;
        }
    }
    
    //optionally correct the high frequencies
    if (option == 1) {
        int ii2, jj2;
        if ( w_test*2 == w ) {
            ii = ind_w; //positive and location in input
            ii2= ind_w + (zoom-1)*w; //negative
            for(int j=0; j<h; j++) {
                jj = (j<ind_h) ? j : j + (zoom-1)*h;
                for (int l = 0; l < 2*pd; l++) {
                    outhat2[jj][ii2][l]  *= 0.5;
                    outhat2[jj][ii][l]  = outhat2[jj][ii2][l];
                }
            }
        }
        
        if ( h_test*2 == h ) {
            jj = ind_h; //positive and location in input
            jj2= ind_h + (zoom-1)*h; //negative
            for(int i=0; i<w; i++) {
                ii = (i<ind_w) ? i : i + (zoom-1)*w;
                for (int l = 0; l < 2*pd; l++) {
                    outhat2[jj2][ii][l]  *= 0.5;
                    outhat2[jj][ii][l]  = outhat2[jj2][ii][l];
                }
            }
        }
        
        if ( h_test*2 == h && w_test*2 == w) {
            double tmp;
            ii = ind_w; //positive location in output and location in input
            ii2= ind_w + (zoom-1)*w; //negative
            jj = ind_h; //positive location in output and location in input
            jj2= ind_h + (zoom-1)*h; //negative
            for (int l = 0; l < 2*pd; l++) {
                tmp = inhat2[jj][ii][l]*norm;
                outhat2[jj][ii][l] = 0.25*tmp;
                outhat2[jj][ii2][l] = 0.25*tmp;
                outhat2[jj2][ii][l] = 0.25*tmp;
                outhat2[jj2][ii2][l] = 0.25*tmp;
            }
        }
    }
  }
  else {
      for (int i=0; i<2*pd*zoom*zoom*w*h; i++)
                outhat[i]=inhat[i];
  }
}

#endif