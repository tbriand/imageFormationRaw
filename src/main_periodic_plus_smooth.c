// Computation of the periodic plus smooth decomposition of an image

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

#include "iio.h"
#include "fft_core.h"
#include "periodic_plus_smooth.h"

int main(int c, char *v[])
{
    if (c < 3) {
            fprintf(stderr,"usage:\n\t%s  in periodic smooth [zoom]\n", *v);
            //                         0  1  2        3      4
            return EXIT_FAILURE;
    }

    char *filename_in = c > 1 ? v[1] : "-";
    char *filename_periodic = c > 2 ? v[2] : "-";
    char *filename_smooth = c > 3 ? v[3] : "-";
    double zoom = c > 4 ? atof(v[4]) : 1.0;

    // read data and convert to double
    int w, h, pd;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);

    // output size of periodic component
    int wout = zoom*w;
    int hout = zoom*h;
    double *periodic = malloc(wout*hout*pd*sizeof(double));
    double *smooth = malloc(w*h*pd*sizeof(double));
    
    // initialize FFTW
    init_fftw();
    
    // compute the periodic plus smooth decomposition
    periodic_plus_smooth_decomposition(periodic, smooth, in, w, h, pd, zoom);
    
    // write output image
    iio_write_image_double_split(filename_periodic, periodic, wout, hout, pd);
    iio_write_image_double_split(filename_smooth, smooth, w, h, pd);
    
    // free memory
    free(in);
    free(periodic);
    free(smooth);
    clean_fftw();
    
    return EXIT_SUCCESS;
} 