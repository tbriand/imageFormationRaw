// Implementation of an affine transformation of an image dynamic
// It maps [min,max] to [0,value] where by default value = 255

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

#include "iio.h"

double min(double *tab, size_t N)
{
  double min=tab[0];
  for(unsigned int i=1;i<N;++i)
  {
     if (tab[i]<min) min=tab[i];
  }
  return(min);
}

double max(double *tab, size_t N)
{
  double max=tab[0];
  for(unsigned int i=1;i<N;++i)
  {
     if (tab[i]>max) max=tab[i];
  }
  return(max);
}

static void affine0255(double *in, double *out, int w, int h, int pd, int type)
{
  double M = max(in,w*h);

  if( type == 0) {
    double m = min(in,w*h);
    for( int i=0;i<w*h*pd;i++){
      out[i] = 255*(in[i]-m)/(M-m);
    }
  }
  else {
    double factor = type*1.0/M;
    for( int i=0;i<w*h*pd;i++){
      out[i] = factor*in[i];
    }
  }
}

int main(int c, char *v[])
{
    if (c < 3) {
            fprintf(stderr,"usage:\n\t%s in out [max]      \n", *v);
            //                         0 1  2   3
            fprintf(stderr,"max: 0 (affine0255) or value (set the max to this value)\n");
            return EXIT_FAILURE;
    }
    char *filename_in = c > 1 ? v[1] : "-";
    char *filename_out = c > 2 ? v[2] : "-";
    int type = c > 3 ? atoi(v[3]) : 0;

    int w, h, pd;
    double *in = iio_read_image_double_vec(filename_in, &w, &h, &pd);
    
    double *out = malloc(w*h*pd*sizeof*out);
    
    affine0255(in, out, w, h, pd, type);
    
    iio_write_image_double_vec(filename_out, out, w, h, pd);
    
    free(in);
    free(out);
    
    return EXIT_SUCCESS;
}