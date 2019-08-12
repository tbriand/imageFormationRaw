// Split a mosaicked image into four grayscale images according to the RGGB patern
// Not possible to write in PNG format if double image
// Raw images are quantified and can be handled with float format

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

int main(int c, char *v[]) {
    if (c < 3) {
        fprintf(stderr,"usage:\n\t%s in path_out ext\n", *v);
        //                         0 1  2        3
        return EXIT_FAILURE;
    }

    char *filename_in = v[1];
    char *path_out = v[2];
    char *ext = c > 3 ? v[3] : "tiff";

    int w, h, pd;
    float *in = iio_read_image_float_vec(filename_in, &w, &h, &pd);

    if ( pd != 1 ) {
        fprintf(stderr,"Images should have one channel\n");
        return EXIT_FAILURE;
    }

    int w1 = (w+1)/2;
    int h1 = (h+1)/2;
    int w2 = w/2;
    int h2 = h/2;

    // channel 1 (RED)
    uint16_t *out1 = malloc(w1*h1*sizeof*out1);
    for(int j=0; j<h1; j++)
        for(int i=0; i<w1; i++)
            out1[i+w1*j] = in[2*i+w*2*j];

    // channel 2 (GREEN)
    uint16_t *out2 = malloc(w2*h1*sizeof*out2);
    for(int j=0; j<h1; j++)
        for(int i=0; i<w2; i++)
            out2[i+w2*j] = in[2*i+1+w*2*j];
    
    // channel 3 (GREEN)
    uint16_t *out3 = malloc(w1*h2*sizeof*out3);
    for(int j=0; j<h2; j++)
        for(int i=0; i<w1; i++)
            out3[i+w1*j] = in[2*i+w*(2*j+1)];
        
    // channel 4 (BLUE)
    uint16_t *out4 = malloc(w2*h2*sizeof*out4);
    for(int j=0; j<h2; j++)
        for(int i=0; i<w2; i++)
            out4[i+w2*j] = in[2*i+1+w*(2*j+1)];

    char filename_out[500];

    sprintf(filename_out,"%s_%i.%s", path_out, 1, ext);
    iio_write_image_uint16_vec(filename_out, out1, w1, h1, 1);
    
    sprintf(filename_out,"%s_%i.%s", path_out, 2, ext);
    iio_write_image_uint16_vec(filename_out, out2, w2, h1, 1);

    sprintf(filename_out,"%s_%i.%s", path_out, 3, ext);
    iio_write_image_uint16_vec(filename_out, out3, w1, h2, 1);    

    sprintf(filename_out,"%s_%i.%s", path_out, 4, ext);
    iio_write_image_uint16_vec(filename_out, out4, w2, h2, 1);
    
    free(in);
    free(out1);
    free(out2);
    free(out3);
    free(out4);

    return EXIT_SUCCESS;
}
