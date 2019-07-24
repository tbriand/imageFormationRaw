#ifndef INTERPOLATION_CORE
#define INTERPOLATION_CORE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "homography_core.c"
#include "pComponent_core.c"
#include "splinter.c"
#include "bicubic_interpolation.c"

#ifndef BOUNDARY_DEFINITION
#define BOUNDARY_DEFINITION
typedef enum
{
    BOUNDARY_CONSTANT = 0,
    BOUNDARY_HSYMMETRIC = 1,
    BOUNDARY_WSYMMETRIC = 2,
    BOUNDARY_PERIODIC = 3
} BoundaryExt;
#endif

// Read boundary extension
static BoundaryExt read_ext(const char* boundary) {
    if(0 == strncmp(boundary, "constant", strlen(boundary)))
        return BOUNDARY_CONSTANT;
    if(0 == strncmp(boundary, "per", strlen(boundary)))
        return BOUNDARY_PERIODIC;
    if(0 == strncmp(boundary, "hsymmetric", strlen(boundary)))
        return BOUNDARY_HSYMMETRIC;
    if(0 == strncmp(boundary, "wsymmetric", strlen(boundary)))
        return BOUNDARY_WSYMMETRIC;
    fprintf(stderr,"Unknown boundary condition %s\n",boundary);
    exit(EXIT_FAILURE);
}

void splinter_at(double *out, double *in, int w, int h, int pd, int n, BoundaryExt bc, 
                 double precision, int larger, double *x, double *y, int numPixels) {
    //init plan (prefiltering)
    splinter_plan_t plan = splinter_plan(in, w, h, pd, n, bc, precision, larger);
    
    //computation of the pixel locations
    double *outp = malloc(pd*sizeof*outp);
    for(int i=0; i<numPixels; i++) {
            splinter(outp, x[i], y[i], plan);
            for(int k=0; k<pd; k++)
                out[k*numPixels] = outp[k];
            ++out;
    }
    
    free(outp);
    splinter_destroy_plan(plan);
}

void interpolate_at(double *out, double *in, int w, int h, int pd, char *interp, BoundaryExt bc, double *x, double *y, int numPixels) {

    if (0 == strncmp(interp, "bicubic", strlen(interp))) {
        for( int l = 0; l < pd; l++ )
            interpolate_bicubic(out + l*numPixels, in + l*w*h, w, h, bc, x, y, numPixels);
    }
    else if ( 0 == strcmp(interp, "spline") ) {
        // eps <-> 10^(-eps)
        double precision = 1;
        for (int i=0; i<12; i++)
            precision *= 0.1;
        // larger algorithm or not
        int larger = 0;
        if ( bc == BOUNDARY_CONSTANT )
            larger = 1;
        splinter_at(out, in, w, h, pd, 11,
                       bc, precision, larger, x, y, numPixels);
    }
    else
        fprintf(stderr,"Unknown interpolation method...\n");
}

void interpolate_at_vec(double *out, double *in, int w, int h, int pd, char *interp, BoundaryExt bc, double *x, double *y, int numPixels) {
    // storage of the input data is not compatible with the interpolation structure
    // in/out in vec
    // buffs in split
    double *buff = malloc(w*h*pd*sizeof*buff);
    double *buff_out = malloc(numPixels*pd*sizeof*buff_out);
    
    // goes from vec to split
    for(int l=0; l<pd; l++)
        for(int i=0; i<w*h; i++)
            buff[i+l*w*h] = in[i*pd+l];
    
    interpolate_at(buff_out, buff, w, h, pd, interp, bc, x, y, numPixels);
    free(buff);
    
    //goes from split to vec
    for(int l=0; l<pd; l++)
        for (int i=0; i<numPixels; i++) 
            out[i*pd+l] = buff_out[i+l*numPixels];
        
    free(buff_out);
}

static void resampling(double *out, double *in, int w, int h, int pd, double H[9], char *interp, BoundaryExt bc, float zoom) {
    
    // output sizes
    int wout = w/zoom;
    int hout = h/zoom;
    
    // create pixel locations
    double iH[9];
    invert_homography(iH, H);
    double *x = malloc(wout*hout*sizeof*x);
    double *y = malloc(wout*hout*sizeof*y);
    double p[2], q[2];
    for (int j = 0; j < hout; j++) {
        p[1] = j*zoom;
        for (int i = 0; i < wout; i++) {
            p[0] = i*zoom;
            apply_homo(q, p, iH);
            x[j*wout+i] = q[0];
            y[j*wout+i] = q[1];
        }
    }
    
    if (0 == strcmp(interp, "per") || 0 == strcmp(interp, "splineper"  )) {
        char *interpmethod;
        if( 0 == strcmp(interp, "splineper" ) ) interpmethod = "spline";
        if( 0 == strcmp(interp, "per" ) )       interpmethod = "bicubic";
        
        int zoomper = 2;
        int wper = w*zoomper;
        int hper = h*zoomper;
        
        // periodic plus smooth decomposition
        double *pComp_in_zoomed = malloc(wper*hper*pd*sizeof(double));
        double *sComp_in = malloc(w*h*pd*sizeof(double));
        compute_component(in, pComp_in_zoomed, sComp_in, zoomper, w, h, pd);
        
        // interpolate smooth
        interpolate_at_vec(out, sComp_in, w, h, pd, interpmethod, bc, x, y, wout*hout);
        free(sComp_in);
        
        // create pixel locations for the zoomed periodic version
        for (int i = 0; i < wout*hout; i++) {
                x[i] *= zoomper;
                y[i] *= zoomper;
        }
        
        // interpolate periodic component
        double *pComp = malloc(wout*hout*pd*sizeof(double));
        interpolate_at_vec(pComp, pComp_in_zoomed, wper, hper, pd, interpmethod, BOUNDARY_PERIODIC, x, y, wout*hout);
        free(pComp_in_zoomed);
        
        // sum
        for (int k = 0; k < wout*hout*pd; k++)
            out[k] += pComp[k];
        free(pComp);
    }
    else {
        interpolate_at_vec(out, in, w, h, pd, interp, bc, x, y, wout*hout);
    }
    
    // free memory
    free(x);
    free(y);
}

#endif