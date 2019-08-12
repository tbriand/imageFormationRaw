// Functions for the equalizing the dynamics of two images

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

#ifndef EQUALIZATION_CORE
#define EQUALIZATION_CORE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* histogram equalization part */

static int f_compare (const void * a, const void * b)
{
    if(*(const double*)a < *(const double*)b)
        return -1;
    return *(const double*)a > *(const double*)b;
}

static void equalization_histo(double *ref, double *modified, double *out, int w, int h, int pd)
{
    // mix values and indexes to keep track of pixels' location
    double *sort_values_ref=malloc(2*w*h*sizeof*sort_values_ref);
    double *sort_values_modified=malloc(2*w*h*sizeof*sort_values_modified);

    for (int l=0; l<pd; l++) {
        for (int idx=0; idx<w*h; idx++) {
            sort_values_ref[2*idx] = ref[idx + l*w*h];
            sort_values_ref[2*idx+1] = (double) idx;
            sort_values_modified[2*idx] = modified[idx + l*w*h];
            sort_values_modified[2*idx+1] = (double) idx;
        }

        // sort pixels depending on their values
        qsort(sort_values_ref, w*h, 2*sizeof(double), f_compare);
        qsort(sort_values_modified, w*h, 2*sizeof(double), f_compare);

        // histogram matching
        for(int idx = 0; idx < w*h ; idx++)
            out[ (int) sort_values_modified[2*idx+1] + l*w*h] = sort_values_ref[2*idx];
    }

    // free memory
    free(sort_values_ref);
    free(sort_values_modified);
}

/* end histogram equalization part */

static void equalization_meanp(double *ref, double *modified, double *out,
                               int w, int h, int w2, int h2, int pd, int raw)
{
    if ( raw == 1 ) {
        double mean_ref[3] = {0};
        double mean_modified[3]={0};
        int count[3] = {0};
        int count2[3] = {0};
        int l;
        
        for(int j = 0; j < h; j++)
            for(int i=0;i<w;i++) {
                l = i%2+j%2;
                mean_ref[l] += ref[i + j*w];
                count[l]++;
            }
        
        for( int j = 0; j < h2; j++)
            for(int i = 0; i < w2; i++) {
                l = i%2+j%2;
                mean_modified[l] += modified[i + j*w2];
                count2[l]++;
            }

        mean_ref[0] /= (double) count[0];
        mean_ref[1] /= (double) count[1];
        mean_ref[2] /= (double) count[2];

        mean_modified[0] /= (double) count2[0];
        mean_modified[1] /= (double) count2[1];
        mean_modified[2] /= (double) count2[2];

        for( int j = 0; j < h2; j++)
            for(int i = 0; i < w2; i++) {
                l=i%2+j%2;
                out[i+j*w2] = modified[i + j*w2] - mean_modified[l] + mean_ref[l];
            }
    }
    else {
        double *mean_ref = malloc(pd*sizeof(double));
        double *mean_modified = malloc(pd*sizeof(double));

        // mean computation
        for(int l = 0; l < pd; l++) {
            mean_ref[l] = 0;
            for(int i = 0; i<w*h; i++)
                mean_ref[l] += ref[i + l*w*h];
            mean_ref[l] /= (double) (w*h);

            mean_modified[l]=0;
            for(int i = 0; i < w2*h2; i++)
                mean_modified[l] += modified[i+l*w2*h2];
            mean_modified[l] /= (double) (w2*h2);
        }

        // equalization by shifting the mean
        for(int l = 0; l < pd; l++)
            for(int i = 0; i < w2*h2; i++)
                out[i + l*w2*h2] = modified[i + l*w2*h2] - mean_modified[l] + mean_ref[l];

        free(mean_ref);
        free(mean_modified);
    }
}

static void equalization_meanx(double *ref, double *modified, double *out,
                                    int w, int h, int w2, int h2, int pd, int raw)
{
    if ( raw == 1 ) {
        double mean_ref[3] = {0};
        double mean_modified[3]={0};
        int count[3] = {0};
        int count2[3] = {0};
        int l;
        
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                l = i%2+j%2;
                mean_ref[l] += ref[i+j*w];
                count[l]++;
            }
        for( int j=0; j<h2; j++)  
            for(int i=0; i<w2; i++) {
                l = i%2+j%2;
                mean_modified[l] += modified[i+j*w];
                count2[l]++;
            }

        mean_ref[0] /= (double) count[0];
        mean_ref[1] /= (double) count[1];
        mean_ref[2] /= (double) count[2];

        mean_modified[0] /= (double) count2[0];
        mean_modified[1] /= (double) count2[1];
        mean_modified[2] /= (double) count2[2];


        double factor[3];
        for(int l=0; l<3; l++)
            factor[l] = mean_ref[l]/mean_modified[l];

        for( int j=0; j<h2; j++)
            for(int i=0; i<w2; i++) {
                    l=i%2+j%2;
                    out[i+j*w2] = modified[i+j*w2]*factor[l];
            }
    }
    else {
        double *mean_ref = malloc(pd*sizeof(double));
        double *mean_modified = malloc(pd*sizeof(double));

        /* mean computation */
        for(int l=0; l<pd; l++) {
            mean_ref[l] = 0;
            for(int i=0; i<w*h; i++)
                    mean_ref[l] += ref[i+l*w*h];
            mean_ref[l] /= (double) (w*h);

            mean_modified[l] = 0;
            for(int i=0; i<w2*h2; i++)
                mean_modified[l] += modified[i+l*w2*h2];
            mean_modified[l] /= (double) (w2*h2);
        }

        /* equalization by shifting the mean */
        for(int l = 0; l < pd; l++) {
            double factor = mean_ref[l]/mean_modified[l];
            for(int i=0; i<w2*h2; i++)
                out[i+l*w2*h2] = modified[i+l*w2*h2]*factor;
        }

        free(mean_ref);
        free(mean_modified);
    }
}

static void equalization_meanx_double(double *ref, double *modified, double *out,
                                    int w, int h, int w2, int h2, int pd, int raw)
{
    if ( raw == 1 ) {
        double mean_ref[3] = {0};
        double mean_modified[3]={0};
        int count[3] = {0};
        int count2[3] = {0};
        int l;
        
        for( int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                l = i%2+j%2;
                mean_ref[l] += ref[i+j*w];
                count[l]++;
            }
        
        for( int j=0; j<h2; j++) 
            for(int i=0; i<w2; i++) {
                l = i%2+j%2;
                mean_modified[l] += modified[i+j*w2];
                count2[l]++;
            }

        mean_ref[0] /= (double) count[0];
        mean_ref[1] /= (double) count[1];
        mean_ref[2] /= (double) count[2];

        mean_modified[0] /= (double) count2[0];
        mean_modified[1] /= (double) count2[1];
        mean_modified[2] /= (double) count2[2];


        double factor[3];
        for(int l=0; l<3; l++)
            factor[l] = mean_ref[l]/mean_modified[l];
        for(int j=0; j<h2; j++)
            for(int i=0; i<w2; i++)
            {
                l=i%2+j%2;
                out[i+j*w2] = modified[i+j*w2+l*w2*h2]*factor[l];
            }
        }
    else {
        double *mean_ref = malloc(pd*sizeof(double));
        double *mean_modified = malloc(pd*sizeof(double));

        /* mean computation */
        for(int l=0; l<pd; l++) {
            mean_ref[l] = 0;
            for(int i=0; i<w*h; i++)
                mean_ref[l] += ref[i+l*w*h];
            mean_ref[l] /= (double) (w*h);

            mean_modified[l]=0;
            for(int i=0; i<w2*h2; i++)
                mean_modified[l] += modified[i+l*w2*h2];
            mean_modified[l] /= (double) (w2*h2);
        }

        /* equalization by shifting the mean */
        for(int l=0; l<pd; l++) {
            double factor = mean_ref[l]/mean_modified[l];
            for(int i=0; i<w2*h2; i++)
                out[i+l*w2*h2] = modified[i+l*w2*h2]*factor;
        }
        
        // free memory
        free(mean_ref);
        free(mean_modified);
    }
}

static void equalization_affine(double *ref, double *modified, double *out,
                                    int w, int h, int w2, int h2, int pd, int raw)
{
    double tmp;

    if ( raw == 1 ) {
        double mean_ref[3] = {0};
        double mean_modified[3]={0};
        double v_ref[3] = {0};
        double v_modified[3] = {0};
        int count[3] = {0};
        int count2[3] = {0};
        int l;
        
        for( int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                l = i%2+j%2;
                mean_ref[l] += ref[i+j*w];
                count[l]++;
            }
        for( int j=0; j<h2; j++)
            for(int i=0;i<w2;i++) {
                l = i%2+j%2;
                mean_modified[l] += modified[i+j*w2];
                count2[l]++;
            }

        mean_ref[0] /= (double) count[0];
        mean_ref[1] /= (double) count[1];
        mean_ref[2] /= (double) count[2];

        mean_modified[0] /= (double) count2[0];
        mean_modified[1] /= (double) count2[1];
        mean_modified[2] /= (double) count2[2];

        for( int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                l = i%2+j%2;
                tmp = ref[i+j*w] - mean_ref[l];
                v_ref[l] += tmp*tmp;
            }
        for(int j=0; j<h2; j++)
            for(int i=0;i<w2;i++) {
                l = i%2+j%2;
                tmp = modified[i+j*w2] - mean_modified[l];
                v_modified[l] += tmp*tmp;
            }

        double pente[3];
        for (int i=0; i<3; i++)
            pente[i] = sqrt(v_ref[i]/v_modified[i]*count2[i]*1.0/((double) count[i]));

        for(int j=0; j<h2; j++)
            for(int i=0; i<w2; i++) {
                l=i%2+j%2;
                out[i+j*w2] = pente[l]*(modified[i+j*w2] - mean_modified[l]) + mean_ref[l];
            }
    }
    else {
        double *mean_ref = malloc(pd*sizeof(double));
        double *mean_modified = malloc(pd*sizeof(double));
        double *v_ref = malloc(pd*sizeof(double));
        double *v_modified = malloc(pd*sizeof(double));
        double *slope = malloc(pd*sizeof(double));

        /* mean computation */
        for(int l=0; l<pd; l++) {
            mean_ref[l] = 0;
            for(int i=0; i<w*h; i++)
                mean_ref[l] += ref[i+l*w*h];
            mean_ref[l] /= (double) (w*h);

            mean_modified[l]=0;
            for(int i=0; i<w2*h2; i++)
                mean_modified[l] += modified[i+l*w2*h2];
            mean_modified[l] /= (double) (w2*h2);
        }

        /* variance computation and affine equalization */
        for(int l=0; l<pd; l++) {
            v_ref[l] = 0;
            for(int i=0; i<w*h; i++) {
                tmp = ref[i+l*w*h] - mean_ref[l];
                v_ref[l] += tmp*tmp;
            }
            v_ref[l] /= (double) (w*h);

            v_modified[l] = 0;
            for(int i=0; i<w2*h2; i++) {
                tmp = modified[i+l*w2*h2] - mean_modified[l];
                v_modified[l] += tmp*tmp;
            }
            v_modified[l] /= (double) (w2*h2);
            
            slope[l] = sqrt(v_ref[l]/v_modified[l]);

            /* affine equalization */
            for(int i=0; i<w2*h2; i++)
                out[i+l*w2*h2] = slope[l]*(modified[i+l*w2*h2] - mean_modified[l]) + mean_ref[l];
        }

        free(mean_ref);
        free(mean_modified);
        free(v_ref);
        free(v_modified);
        free(slope);
    }
}

void channel_equalization(double *in, double *out, int w, int h, int pd)
{
    if ( pd == 1 ) { //raw case
        double sum_ref = 0.0;
        double sum_in[3]={0.0};
        int count[3] = {0.0};
        int l;
        
        for(int i=0; i<w*h; i++)
            sum_ref += in[i];
        
        for(int j=0; j<h; j++)
            for(int i=0;i<w;i++) {
                l = i%2+j%2;
                sum_in[l] += in[i+j*w];
                count[l]++;
        }

        double factor[3];
        for(int i=0; i<3; i++)
            factor[i] = sum_ref/sum_in[i]*count[i]*1.0/((double) w*h); // normalize by the ration of elements in each channel
        
        for(int j=0; j<h; j++)
            for(int i=0; i<w; i++) {
                l=i%2+j%2;
                out[i+j*w] = in[i+j*w]*factor[l];
            }
    }
    else if ( pd == 3 ) {
        double sum_in[3] = {0.0};

        for (int l=0; l<3; l++)
            for(int i=0; i<w*h; i++)
                sum_in[l] += in[i+l*w*h];

        double sum = 0.0;
        for (int l=0; l<3; l++)
            sum += sum_in[l];
        
        double factor[3];
        for(int l=0; l<3; l++)
            factor[l] = sum/(3*sum_in[l]);
        
        for (int l=0; l<3; l++)
            for(int i=0; i<w*h; i++)
                out[i+l*w*h] = in[i+l*w*h]*factor[l];
    
    }
}

#endif