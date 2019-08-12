// Application of a lowpass filter

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

#include "fft_core.h"
#include "filter_core.h"
#include "iio.h"

int main(int c, char *v[])
{
    if (c != 4) {
            fprintf(stderr,"usage:\n\t%s  in out filter_type\n", *v);
            //                         0  1  2   3
            fprintf(stderr,"filter type : 0 --> perfect lowpass\n");
            fprintf(stderr,"filter type : 1 --> hvs\n");
            fprintf(stderr,"filter type : 2 --> bilinear+luminance\n");
            return EXIT_FAILURE;
    }

    char *filename_in = c > 1 ? v[1] : "-";
    char *filename_out = c > 2 ? v[2] : "-";
    int filter_type = atoi(v[3]);
    
    if ( filter_type < 0 || filter_type > 2 )
        filter_type = 0;

    // read data and convert to double
    int w, h, pd;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
    double *out = malloc(w*h*pd*sizeof*out);

    // initialize FFTW
    init_fftw();
    
    // apply the lowpass filter
    if ( filter_type == 0) // available for every image
        apply_perfect_lowpass(in, out, w, h, pd);
    else {
        if( pd != 1) { // only for mosaicked images
            printf("This filter can only be applied to mosaicked images\n");
            return EXIT_FAILURE;
        }

        if ( filter_type == 1 )
            apply_hvs_lowpass(in, out, w, h);
        else if ( filter_type == 2 )
            apply_bilinear_luminance( in, out, w, h );
    }
    
    // write output image
    iio_write_image_double_split(filename_out, out, w, h, pd);
    
    // free memory
    free(in);
    free(out);
    clean_fftw();
    
    return EXIT_SUCCESS;
} 
