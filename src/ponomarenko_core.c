// utility functions regarding the VST from ponomarenko estimation

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

#ifndef PONOMARENKO_CORE
#define PONOMARENKO_CORE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parsenumbers.c"
#include "xfopen.c"
#include "iio.h"

// regression y = sqrt(A+Bx) i.e y^2 = A + Bx
static void regression_ponomarenko(const double *t, int n, double A[2]) {
    // sanity check
    int nbins = n/2;
    if ( nbins*2 != n ) {
        fprintf(stderr,"Wrong content in the input file\n");
    }
    else {
        // Matlab equivalent
        // X = [ones(nbins,1),)
        // x = t[0:2:n-1];
        // y = t[1:2:n-1]^2;
        // [A;B] = (X'*X)^{-1} X*y;
        // X'*X = [nbins,sum(x);sum(x),sum(x.^2)]
        // X*y = [sum(y);sum(x.*y)]
        double sx = 0;
        double sx2 = 0;
        double sy = 0;
        double sxy = 0;
        for( int i=0; i<nbins; i++) {
            sx += t[2*i];
            sx2 += t[2*i]*t[2*i];
            sy += t[2*i+1]*t[2*i+1];
            sxy += t[2*i]*t[2*i+1]*t[2*i+1];
        }
        
        double idet = nbins*sx2 - sx*sx;
        if ( idet == 0 )
            fprintf(stderr,"Not inversible matrix\n");
        else {
            idet = 1.0/idet;

            A[0] = idet*(sx2*sy - sx*sxy);
            A[1] = idet*(-sx*sy + nbins*sxy);
        
           // printf("%lf %lf\n",A[0],A[1]);
        }
    }
}

static void regression_1channel(char *filename, double C[2])
{
    FILE *f = xfopen(filename, "r");
    int n;
    double *t = read_ascii_doubles(f, &n);
    xfclose(f);
    
    // regression
    regression_ponomarenko(t, n, C);
    free(t);
}

static void regression_3channels(char *inpat, int ind, double R[2], double G[2], double B[2])
{
    char filename[1000];
    
    /* red channel */
    // read file
    if ( ind == -1 )
        sprintf(filename,"%s_1.txt", inpat);
    else 
        sprintf(filename,"%s_%i_1.txt", inpat, ind);
    FILE *fr = xfopen(filename, "r");
    int nr;
    double *r = read_ascii_doubles(fr, &nr);
    xfclose(fr);
    
    // regression
    regression_ponomarenko(r, nr, R);
    free(r);
    
    /* blue channel */
    // read file
    if ( ind == -1 )
        sprintf(filename,"%s_4.txt", inpat);
    else 
        sprintf(filename,"%s_%i_4.txt", inpat, ind);
    FILE *fb = xfopen(filename, "r");
    int nb;
    double *b = read_ascii_doubles(fb, &nb);
    xfclose(fb);
    
    // regression
    regression_ponomarenko(b, nb, B);
    free(b);
    
    /* green channel */
    // read file
    if ( ind == -1 )
        sprintf(filename,"%s_2.txt", inpat);
    else 
        sprintf(filename,"%s_%i_2.txt", inpat, ind);
    FILE *fg1 = xfopen(filename, "r");
    int ng1;
    double *g1 = read_ascii_doubles(fg1, &ng1);
    xfclose(fg1);
    
    if ( ind == -1 )
        sprintf(filename,"%s_3.txt", inpat);
    else 
        sprintf(filename,"%s_%i_3.txt", inpat, ind);
    FILE *fg2 = xfopen(filename, "r");
    int ng2;
    double *g2 = read_ascii_doubles(fg2, &ng2);
    xfclose(fg2);
    
    int ng = ng1 + ng2;
    double *g = malloc(ng*sizeof*g);
    for(int i = 0; i < ng1; i++)
        g[i] = g1[i];
    for(int i = 0; i < ng2; i++)
        g[ng1 + i] = g2[i];
    
    // regression
    regression_ponomarenko(g, ng, G);
    free(g1);
    free(g2);
    free(g);
}

static void vst_1channel(char *filename_in, char *filename_out, const double C[2])
{
    int w, h, pd;
    double *in = iio_read_image_double_vec(filename_in, &w, &h, &pd);
    if ( pd != 1 )
        fprintf(stderr,"Wrong number of channels (pd = %i) for image %s\n", pd, filename_in);
    double *out = malloc( w * h * sizeof*out );

    for( int i = 0; i < w*h; i++) 
        out[i] = sqrt(C[0]+C[1]*in[i]);
    iio_write_image_double_vec(filename_out, out, w, h, 1);
    
    free(in);
    free(out);
}

static void vst_3channels(char *inpat, char *ext, char *outpat, int ind, double R[2], double G[2], double B[2])
{
    char filename_in[1000];
    char filename_out[1000];
 
    if ( ind == -1 ) {
        sprintf(filename_in, "%s", inpat);
        sprintf(filename_out, "%s", outpat);
    }
    else {
        sprintf(filename_in,"%s_%i.%s", inpat, ind, ext);
        sprintf(filename_out,"%s_%i.tiff", outpat, ind);
    }
    
    int w, h, pd;
    double *in = iio_read_image_double_vec(filename_in, &w, &h, &pd);
    if ( pd != 1 )
        fprintf(stderr,"Wrong number of channels (pd = %i) for image %s\n", pd, filename_in);
    double *out = malloc( w * h * sizeof*out );

    int l;
    double A, C;
    for (int j = 0; j < h; j++) {
        for( int i = 0; i < w; i++) {
            l = i%2 + j%2;
            if ( l == 0 ) { // red
                A = R[0];
                C = R[1];
            }
            else if (l == 1 ) { // green
                A = G[0];
                C = G[1];
            }
            else { // blue
                A = B[0];
                C = B[1];
            }
            out[i+w*j] = (A+C*in[i+w*j] >=0 ) ? sqrt(A+C*in[i+w*j]) : 0;
        }
    }
    iio_write_image_double_vec(filename_out, out, w, h, 1);
    
    free(in);
    free(out);
}

static void trapeze_integral_tools(double *t, double *pos, double *val, double *slope, double *integral, int n) {
    // extract position and set first to 0
    pos[0] = 0;
    for (int i = 1; i < n; i++)
        pos[i] = t[2*i-2];
    
    // compute value of the integrand (lack the first value that is obtained thanks to extrapolation later)
    for (int i = 1; i < n; i++)
        val[i] = 1/t[2*i-1];
    
    // compute the slopes (lack the first and last value)
    for (int i = 1; i < n - 1; i++) 
        slope[i] = (val[i+1] - val[i])/(pos[i+1] - pos[i]);
    
    // first and last slopes
    slope[0] = slope[1];
    slope[n-1] = slope[n-2];
    
    // first integrand value obtained by extrapolation but it should be non-negative
    val[0] = val[1] - pos[1]*slope[1];
    if ( val[0] < 0 ) {
        pos[0] = pos[1] - val[1]/slope[1];
        val[0] = 0;
    }
    
    integral[0] = 0.5*(val[1]+val[0])*(pos[1] - pos[0]);
    // compute the integral values as cumulative sums
    for (int i = 1; i < n - 1; i++)
        integral[i] = integral[i-1] + 0.5*(val[i+1]+val[i])*(pos[i+1] - pos[i]);
}
 
// note that the output is necessarily non-negative !
// no need to check even for i =0
static double trapeze_integral(double x, double *pos, double *val, double *slope, double *integral, int n) {
    // find the interval i such that x_i < x < x_{i+1}
    int i = 0;
    if ( x <= pos[i]) {
        return 0;
    }
    else {
        while ( (i < n - 1) && (x > pos[ i + 1]) )
            i++;
        
        // integral of the previous intervals
        double y = ( i > 0 ) ? integral[i-1] : 0;
        // add the remaining part
        double y2 = (x - pos[i])*(val[i] + 0.5*slope[i]*(x - pos[i]));
        //printf("y = %lf et y2 = %lf\n",y,y2);
        y = (y2 > 0 ) ? y + y2 : y; // sanity check but y2 should always be positive
        return y;
    }
}

static void vst_trapeze(char *inpat, char *image_in, char *image_out) {
    
    char filename[1000];
    
    //printf("red channel\n");
    /* red channel */
    sprintf(filename,"%s_1.txt", inpat);
    FILE *fr = xfopen(filename, "r");
    int nr;
    double *r = read_ascii_doubles(fr, &nr);
    xfclose(fr);
    int nr2 = nr/2 + 1; // nr/2 points + the zero
    double *rpos = malloc(nr2*sizeof*rpos);
    double *rval = malloc(nr2*sizeof*rval);
    double *rslope = malloc(nr2*sizeof*rslope);
    double *rint = malloc((nr2-1)*sizeof*rint);
    trapeze_integral_tools(r, rpos, rval, rslope, rint, nr2);
    
    /* blue channel */
    sprintf(filename,"%s_4.txt", inpat);
    FILE *fb = xfopen(filename, "r");
    int nb;
    double *b = read_ascii_doubles(fb, &nb);
    xfclose(fb);
    int nb2 = nb/2 + 1;
    double *bpos = malloc(nb2*sizeof*bpos);
    double *bval = malloc(nb2*sizeof*bval);
    double *bslope = malloc(nb2*sizeof*bslope);
    double *bint = malloc((nb2-1)*sizeof*bint);
    trapeze_integral_tools(b, bpos, bval, bslope, bint, nb2);

    /* green channel */
    sprintf(filename,"%s_2.txt", inpat);
    FILE *fg1 = xfopen(filename, "r");
    int ng1;
    double *g1 = read_ascii_doubles(fg1, &ng1);
    xfclose(fg1);
    
    sprintf(filename,"%s_3.txt", inpat);
    FILE *fg2 = xfopen(filename, "r");
    int ng2;
    double *g2 = read_ascii_doubles(fg2, &ng2);
    xfclose(fg2);
    
    // merge the two arrays but keep it sorted with respect to the position
        int ng = ng1 + ng2;
        double *g = malloc(ng*sizeof*g);
        int ig1 = 0, ig2 = 0, ig = 0;
        while (ig1 < ng1 && ig2 < ng2) {
            if ( g1[ig1] < g2[ig2] ) {
                g[ig++] = g1[ig1++]; // position
                g[ig++] = g1[ig1++]; // value
            }
            else {
                g[ig++] = g2[ig2++]; // position
                g[ig++] = g2[ig2++]; // value
            }
        }
        // one of the two is not in the array yet
        while( ig1 < ng1 ) {
                g[ig++] = g1[ig1++]; // position
                g[ig++] = g1[ig1++]; // value
        }
        while( ig2 < ng2 ) {
                g[ig++] = g2[ig2++]; // position
                g[ig++] = g2[ig2++]; // value
        }

    int ng3 = ng/2 + 1;
    double *gpos = malloc(ng3*sizeof*gpos);
    double *gval = malloc(ng3*sizeof*gval);
    double *gslope = malloc(ng3*sizeof*gslope);
    double *gint = malloc((ng3-1)*sizeof*gint);
    trapeze_integral_tools(g, gpos, gval, gslope, gint, ng3);
    
    // read image
    int w, h, pd;
    double *in = iio_read_image_double_vec(image_in, &w, &h, &pd);
    if ( pd != 1 )
        fprintf(stderr,"Wrong number of channels (pd = %i) for image %s\n", pd, image_in);

    // vst transform of the input
    double *out = malloc( w * h * sizeof*out );
    int l;
    for (int j = 0; j < h; j++) {
        for( int i = 0; i < w; i++) {
            l = i%2 + j%2;
            if ( l == 0 ) // red
                out[i+w*j] = trapeze_integral(in[i+w*j], rpos, rval, rslope, rint, nr2);
            else if (l == 1 ) // green
                out[i+w*j] = trapeze_integral(in[i+w*j], gpos, gval, gslope, gint, ng3);
            else // blue
                out[i+w*j] = trapeze_integral(in[i+w*j], bpos, bval, bslope, bint, nb2);
        }
    }
    
    //printf("Write output\n");
    iio_write_image_double_vec(image_out, out, w, h, 1);
    
    // free memory
    free(r);
    free(rpos);
    free(rval);
    free(rslope);
    free(rint);
    free(b);
    free(bpos);
    free(bval);
    free(bslope);
    free(bint);
    free(g);
    free(gpos);
    free(gval);
    free(gslope);
    free(gint);
    free(g1);
    free(g2);
    free(in);
    free(out);
}

static void vst_trapeze_multiple(char *inpat, char *inpat_img, char *ext, char *outpat_img, int ini, int end) {
    
    char filename[1000];
    
    //printf("red channel\n");
    /* red channel */
    sprintf(filename,"%s_1.txt", inpat);
    FILE *fr = xfopen(filename, "r");
    int nr;
    double *r = read_ascii_doubles(fr, &nr);
    xfclose(fr);
    int nr2 = nr/2 + 1; // nr/2 points + the zero
    double *rpos = malloc(nr2*sizeof*rpos);
    double *rval = malloc(nr2*sizeof*rval);
    double *rslope = malloc(nr2*sizeof*rslope);
    double *rint = malloc((nr2-1)*sizeof*rint);
    trapeze_integral_tools(r, rpos, rval, rslope, rint, nr2);
    
    /* blue channel */
    sprintf(filename,"%s_4.txt", inpat);
    FILE *fb = xfopen(filename, "r");
    int nb;
    double *b = read_ascii_doubles(fb, &nb);
    xfclose(fb);
    int nb2 = nb/2 + 1;
    double *bpos = malloc(nb2*sizeof*bpos);
    double *bval = malloc(nb2*sizeof*bval);
    double *bslope = malloc(nb2*sizeof*bslope);
    double *bint = malloc((nb2-1)*sizeof*bint);
    trapeze_integral_tools(b, bpos, bval, bslope, bint, nb2);
    
    /* green channel */
    sprintf(filename,"%s_2.txt", inpat);
    FILE *fg1 = xfopen(filename, "r");
    int ng1;
    double *g1 = read_ascii_doubles(fg1, &ng1);
    xfclose(fg1);
    
    sprintf(filename,"%s_3.txt", inpat);
    FILE *fg2 = xfopen(filename, "r");
    int ng2;
    double *g2 = read_ascii_doubles(fg2, &ng2);
    xfclose(fg2);
    
    // merge the two arrays but keep it sorted with respect to the position
        int ng = ng1 + ng2;
        double *g = malloc(ng*sizeof*g);
        int ig1 = 0, ig2 = 0, ig = 0;
        while (ig1 < ng1 && ig2 < ng2) {
            if ( g1[ig1] < g2[ig2] ) {
                g[ig++] = g1[ig1++]; // position
                g[ig++] = g1[ig1++]; // value
            }
            else {
                g[ig++] = g2[ig2++]; // position
                g[ig++] = g2[ig2++]; // value
            }
        }
        // one of the two is not in the array yet
        while( ig1 < ng1 ) {
                g[ig++] = g1[ig1++]; // position
                g[ig++] = g1[ig1++]; // value
        }
        while( ig2 < ng2 ) {
                g[ig++] = g2[ig2++]; // position
                g[ig++] = g2[ig2++]; // value
        }

    int ng3 = ng/2 + 1;
    double *gpos = malloc(ng3*sizeof*gpos);
    double *gval = malloc(ng3*sizeof*gval);
    double *gslope = malloc(ng3*sizeof*gslope);
    double *gint = malloc((ng3-1)*sizeof*gint);
    trapeze_integral_tools(g, gpos, gval, gslope, gint, ng3);
    
    // transform images
    int w, h, pd;
    double *img_in, *img_out;
    char filename_img_in[1000], filename_img_out[1000];
    for (int index = ini; index <= end; index++) {
        sprintf(filename_img_in,"%s_%i.%s", inpat_img, index, ext);
        img_in = iio_read_image_double_vec(filename_img_in, &w, &h, &pd);
        if ( pd != 1 )
            fprintf(stderr,"Wrong number of channels (pd = %i) for image %s\n", pd, filename_img_in);
        
        // vst transform of the input
        img_out = malloc( w * h * sizeof*img_out );
        int l;
        for (int j = 0; j < h; j++) {
            for( int i = 0; i < w; i++) {
                l = i%2 + j%2;
                if ( l == 0 ) // red
                    img_out[i+w*j] = trapeze_integral(img_in[i+w*j], rpos, rval, rslope, rint, nr2);
                else if (l == 1 ) // green
                    img_out[i+w*j] = trapeze_integral(img_in[i+w*j], gpos, gval, gslope, gint, ng3);
                else // blue
                    img_out[i+w*j] = trapeze_integral(img_in[i+w*j], bpos, bval, bslope, bint, nb2);
            }
        }
        
        sprintf(filename_img_out,"%s_%i.tiff", outpat_img, index);
        iio_write_image_double_vec(filename_img_out, img_out, w, h, 1);
        free(img_in);
        free(img_out);
    }
    
    // free memory
    free(r);
    free(rpos);
    free(rval);
    free(rslope);
    free(rint);
    free(b);
    free(bpos);
    free(bval);
    free(bslope);
    free(bint);
    free(g);
    free(gpos);
    free(gval);
    free(gslope);
    free(gint);
    free(g1);
    free(g2);
 }
 
#endif