// Code for applying the asymptotic equivalent filter (or its inverse)
// of the classical kernel regression to an image using DCT.

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
#include <stdlib.h>
#include <stdio.h>
#include "iio.h"
#include "fft_core.h"
#include "filter_core.h"

int main(int c, char *v[])
{
    if (c != 6) {
            fprintf(stderr,"usage:\n\t%s in    out     order sigma inverse\n", *v);
            //                         0 1     2       3     4     5
            fprintf(stderr,"Apply the asymptotic equivalent filter (or its inverse) to an image using DCT\n");
            return EXIT_FAILURE;
    }
    char *filename_in = c > 1 ? v[1] : "-";
    char *filename_out = c > 2 ? v[2] : "-";
    int order = atoi(v[3]);
    double sigma = atof(v[4]);
    int inverse = atoi(v[5]);
    
    //sanity check
    if ( order != 0 && order != 1 && order != 2 ) {
        fprintf(stderr,"Order value must be 0,1 or 2 (here order=%i)", order);
        return EXIT_FAILURE;
    }
    
    if ( inverse != 0 && inverse != 1 )
        inverse = 0;
    
    // read data
    int w, h, pd;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
    double *out = malloc(w*h*pd*sizeof*out);
    
    // initialize FFTW
    init_fftw();
    
    // apply filter
    apply_asymptotic_nc_filter_global(in, out, w, h, pd, order, sigma, inverse);
    
    // write output image
    iio_write_image_double_split(filename_out, out, w, h, pd);
    
    // free memory
    free(in);
    free(out);
    clean_fftw();
    
    return EXIT_SUCCESS;
}
