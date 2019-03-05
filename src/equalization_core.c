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

int f_compare (const void * a, const void * b)
{
    if(*(const double*)a < *(const double*)b)
        return -1;
    return *(const double*)a > *(const double*)b;
}

static void equalization_histo(double *ref, double *modified, double *out, int w, int h, int pd)
{
  double (*ref2)[pd] = (void*)ref;
  double (*modified2)[pd] = (void*)modified;
  double (*out2)[pd] = (void*)out;

  /* mix values and indexes to keep track of pixels' location */
  double *sort_values_ref=malloc(2*w*h*sizeof*sort_values_ref);
  double *sort_values_modified=malloc(2*w*h*sizeof*sort_values_modified);

  for (int l=0; l<pd;l++) {
    for (int idx=0; idx<w*h; idx++) {
      sort_values_ref[2*idx] = ref2[idx][l];
      sort_values_ref[2*idx+1] = (double) idx;
      sort_values_modified[2*idx] = modified2[idx][l];
      sort_values_modified[2*idx+1] = (double) idx;
    }

    /* sort pixels depending on their values*/
    qsort(sort_values_ref, w*h, 2*sizeof(double), f_compare);
    qsort(sort_values_modified, w*h, 2*sizeof(double), f_compare);

    /* histogram matching */
    for(int idx=0; idx < w*h ; idx++)
      out2[ (int) sort_values_modified[2*idx+1]][l]=sort_values_ref[2*idx];
  }

  /* free memory */
  free(sort_values_ref);
  free(sort_values_modified);
}

/* end histogram equalization part */

static void equalization_meanp(double *ref, double *modified, double *out,
                               int w, int h, int w2, int h2, int pd, int raw)
{
  if ( raw == 1 ) {
    double (*ref2)[w] = (void *) ref;
    double (*modified2)[w2] = (void *) modified;
    double (*out2)[w2] = (void *) out;

    double mean_ref[3] = {0};
    double mean_modified[3]={0};
    int count[3] = {0};
    int count2[3] = {0};
    int l;
    
    for(int i=0;i<w;i++)
      for( int j=0; j<h; j++) {
        l = i%2+j%2;
        mean_ref[l] += ref2[j][i];
        count[l]++;
      }
    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l = i%2+j%2;
        mean_modified[l] += modified2[j][i];
        count2[l]++;
      }

    mean_ref[0] /= (double) count[0];
    mean_ref[1] /= (double) count[1];
    mean_ref[2] /= (double) count[2];

    mean_modified[0] /= (double) count2[0];
    mean_modified[1] /= (double) count2[1];
    mean_modified[2] /= (double) count2[2];

    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l=i%2+j%2;
        out2[j][i] = modified2[j][i]-mean_modified[l]+mean_ref[l];
      }
  }
  else {
    double (*ref2)[w][pd] = (void *) ref;
    double (*modified2)[w2][pd] = (void *) modified;
    double (*out2)[w2][pd] = (void *) out;

    double *mean_ref = malloc(pd*sizeof(double));
    double *mean_modified = malloc(pd*sizeof(double));

    /* mean computation */
    for(int l=0; l<pd; l++) {
      mean_ref[l] = 0;
      for(int i=0; i<w; i++)
        for( int j=0; j<h; j++)
          mean_ref[l] += ref2[j][i][l];
      mean_ref[l] /= (double) (w*h);

      mean_modified[l]=0;
      for(int i=0;i<w2;i++)
        for( int j=0; j<h2; j++)
          mean_modified[l] += modified2[j][i][l];
      mean_modified[l] /= (double) (w2*h2);
    }

    /* equalization by shifting the mean */
    for(int l=0; l<pd; l++) {
      for(int i=0;i<w2;i++) {
        for(int j=0; j<h2; j++) {
          out2[j][i][l] = modified2[j][i][l]-mean_modified[l]+mean_ref[l];
        }
      }
    }

    free(mean_ref);
    free(mean_modified);
  }
}

static void equalization_meanx(double *ref, double *modified, double *out,
                                    int w, int h, int w2, int h2, int pd, int raw)
{
  if ( raw == 1 ) {
    double (*ref2)[w] = (void *) ref;
    double (*modified2)[w2] = (void *) modified;
    double (*out2)[w2] = (void *) out;

    double mean_ref[3] = {0};
    double mean_modified[3]={0};
    int count[3] = {0};
    int count2[3] = {0};
    int l;
    
    for(int i=0;i<w;i++)
      for( int j=0; j<h; j++) {
        l = i%2+j%2;
        mean_ref[l] += ref2[j][i];
        count[l]++;
      }
    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l = i%2+j%2;
        mean_modified[l] += modified2[j][i];
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

    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l=i%2+j%2;
        out2[j][i] = modified2[j][i]*factor[l];
      }
  }
  else {
    double (*ref2)[w][pd] = (void *) ref;
    double (*modified2)[w2][pd] = (void *) modified;
    double (*out2)[w2][pd] = (void *) out;

    double *mean_ref = malloc(pd*sizeof(double));
    double *mean_modified = malloc(pd*sizeof(double));

    /* mean computation */
    for(int l=0; l<pd; l++) {
      mean_ref[l] = 0;
      for(int i=0; i<w; i++)
        for( int j=0; j<h; j++)
          mean_ref[l] += ref2[j][i][l];
      mean_ref[l] /= (double) (w*h);

      mean_modified[l]=0;
      for(int i=0;i<w2;i++)
        for( int j=0; j<h2; j++)
          mean_modified[l] += modified2[j][i][l];
      mean_modified[l] /= (double) (w2*h2);
    }

    /* equalization by shifting the mean */
    for(int l=0; l<pd; l++) {
      double factor = mean_ref[l]/mean_modified[l];
      for(int i=0;i<w2;i++) {
        for(int j=0; j<h2; j++) {
          out2[j][i][l] = modified2[j][i][l]*factor;
        }
      }
    }

    free(mean_ref);
    free(mean_modified);
  }
}

static void equalization_meanx_double(double *ref, double *modified, double *out,
                                    int w, int h, int w2, int h2, int pd, int raw)
{
  if ( raw == 1 ) {
    double (*ref2)[w] = (void *) ref;
    double (*modified2)[w2] = (void *) modified;
    double (*out2)[w2] = (void *) out;

    double mean_ref[3] = {0};
    double mean_modified[3]={0};
    int count[3] = {0};
    int count2[3] = {0};
    int l;
    
    for(int i=0;i<w;i++)
      for( int j=0; j<h; j++) {
        l = i%2+j%2;
        mean_ref[l] += ref2[j][i];
        count[l]++;
      }
    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l = i%2+j%2;
        mean_modified[l] += modified2[j][i];
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
    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l=i%2+j%2;
        out2[j][i] = modified2[j][i]*factor[l];
      }
  }
  else {
    double (*ref2)[w][pd] = (void *) ref;
    double (*modified2)[w2][pd] = (void *) modified;
    double (*out2)[w2][pd] = (void *) out;

    double *mean_ref = malloc(pd*sizeof(double));
    double *mean_modified = malloc(pd*sizeof(double));

    /* mean computation */
    for(int l=0; l<pd; l++) {
      mean_ref[l] = 0;
      for(int i=0; i<w; i++)
        for( int j=0; j<h; j++)
          mean_ref[l] += ref2[j][i][l];
      mean_ref[l] /= (double) (w*h);

      mean_modified[l]=0;
      for(int i=0;i<w2;i++)
        for( int j=0; j<h2; j++)
          mean_modified[l] += modified2[j][i][l];
      mean_modified[l] /= (double) (w2*h2);
    }

    /* equalization by shifting the mean */
    for(int l=0; l<pd; l++) {
      double factor = mean_ref[l]/mean_modified[l];
      for(int i=0;i<w2;i++) {
        for(int j=0; j<h2; j++) {
          out2[j][i][l] = modified2[j][i][l]*factor;
        }
      }
    }

    free(mean_ref);
    free(mean_modified);
  }
}

static void equalization_affine(double *ref, double *modified, double *out,
                                    int w, int h, int w2, int h2, int pd, int raw)
{
  if ( raw == 1 ) {
    double (*ref2)[w] = (void *) ref;
    double (*modified2)[w2] = (void *) modified;
    double (*out2)[w2] = (void *) out;

    double mean_ref[3] = {0};
    double mean_modified[3]={0};
    double v_ref[3] = {0};
    double v_modified[3] = {0};
    int count[3] = {0};
    int count2[3] = {0};
    int l;
    
    for(int i=0;i<w;i++)
      for( int j=0; j<h; j++) {
        l = i%2+j%2;
        mean_ref[l] += ref2[j][i];
        count[l]++;
      }
    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l = i%2+j%2;
        mean_modified[l] += modified2[j][i];
        count2[l]++;
      }

    mean_ref[0] /= (double) count[0];
    mean_ref[1] /= (double) count[1];
    mean_ref[2] /= (double) count[2];

    mean_modified[0] /= (double) count2[0];
    mean_modified[1] /= (double) count2[1];
    mean_modified[2] /= (double) count2[2];

    for(int i=0;i<w;i++)
      for( int j=0; j<h; j++)
        v_ref[i%2+j%2] += (ref2[j][i]-mean_ref[i%2+j%2])*(ref2[j][i]-mean_ref[i%2+j%2]);
    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++)
        v_modified[i%2+j%2] += (modified2[j][i]-mean_modified[i%2+j%2])*(modified2[j][i]-mean_modified[i%2+j%2]);

    double pente[3];
    for (int i=0;i<3;i++) pente[i] = sqrt(v_ref[i]/v_modified[i]*count2[i]*1.0/((double) count[i]));

    for(int i=0;i<w2;i++)
      for( int j=0; j<h2; j++) {
        l=i%2+j%2;
        out2[j][i] = pente[l]*(modified2[j][i]-mean_modified[l])+mean_ref[l];
      }
  }
  else {
    double (*ref2)[w][pd] = (void *) ref;
    double (*modified2)[w2][pd] = (void *) modified;
    double (*out2)[w2][pd] = (void *) out;

    double *mean_ref = malloc(pd*sizeof(double));
    double *mean_modified = malloc(pd*sizeof(double));
    double *v_ref = malloc(pd*sizeof(double));
    double *v_modified = malloc(pd*sizeof(double));
    double *slope = malloc(pd*sizeof(double));

    /* mean computation */
    for(int l=0; l<pd; l++) {
      mean_ref[l] = 0;
      for(int i=0; i<w; i++)
        for( int j=0; j<h; j++)
          mean_ref[l] += ref2[j][i][l];
      mean_ref[l] /= (double) (w*h);

      mean_modified[l]=0;
      for(int i=0;i<w2;i++)
        for( int j=0; j<h2; j++)
          mean_modified[l] += modified2[j][i][l];
      mean_modified[l] /= (double) (w2*h2);
    }

    /* variance computation and affine equalization */
    double tmp;
    for(int l=0; l<pd; l++) {
        
      v_ref[l] = 0;
      for(int i=0;i<w;i++)
        for( int j=0; j<h; j++) {
          tmp = ref2[j][i][l]-mean_ref[l];
          v_ref[l] += tmp*tmp;
        }
      v_ref[l] /= (double) (w*h);

      v_modified[l] = 0;
      for(int i=0;i<w2;i++)
        for( int j=0; j<h2; j++) {
          tmp = modified2[j][i][l]-mean_modified[l];
          v_modified[l] += tmp*tmp;
        }
      v_modified[l] /= (double) (w2*h2);
      

      slope[l] = sqrt(v_ref[l]/v_modified[l]);

      /* affine equalization */
      for(int i=0;i<w2;i++) {
        for(int j=0; j<h2; j++) {
          out2[j][i][l] = slope[l]*(modified2[j][i][l]-mean_modified[l])+mean_ref[l];
        }
      }
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
    double (*in2)[w] = (void *) in;
    double (*out2)[w] = (void *) out;

    double sum_ref = 0;
    double sum_in[3]={0};
    int count[3] = {0};
    int l;
    
    for(int i=0;i<w;i++)
        for( int j=0; j<h; j++)
        sum_ref += in2[j][i];
    for(int i=0;i<w;i++)
        for( int j=0; j<h; j++) {
        l = i%2+j%2;
        sum_in[l] += in2[j][i];
        count[l]++;
        }

    double factor[3];
    for(int i=0; i<3; i++)
        factor[i] = sum_ref/sum_in[i]*count[i]*1.0/((double) w*h); // normalize by the ration of elements in each channel
    
    for(int i=0;i<w;i++)
        for( int j=0; j<h; j++) {
        l=i%2+j%2;
        out2[j][i] = in2[j][i]*factor[l];
        }
  }
  else if ( pd == 3 ) {
    double (*in2)[w][3] = (void *) in;
    double (*out2)[w][3] = (void *) out;

    double sum_in[3]={0};

    for(int i=0;i<w;i++)
        for( int j=0; j<h; j++)
            for ( int l=0; l<3; l++)
                sum_in[l] += in2[j][i][l];

    double sum = 0;
    for (int l=0; l<3; l++)
        sum += sum_in[l];
    
    double factor[3];
    for(int l=0; l<3; l++)
        factor[l] = sum/(3*sum_in[l]);
    
    for(int i=0;i<w;i++)
        for( int j=0; j<h; j++)
            for ( int l=0; l<3; l++)
                out2[j][i][l] = in2[j][i][l]*factor[l];
    
  }
}

#endif