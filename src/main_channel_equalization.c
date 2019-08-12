// Multiplicative mean equalization of the channels of an image.
// When the image has one channel it is assumed to be a mosaicked image
// with a RGGB pattern.

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
#include "equalization_core.h"
#include "iio.h"

int main(int c, char *v[])
{
    if (c < 3) {
        fprintf(stderr,"usage:\n\t%s in  out\n", *v);
        //                         0 1   2
        return EXIT_FAILURE;
    }

    char *filename_in = c > 1 ? v[1] : "-";
    char *filename_out = c > 2 ? v[2] : "-";

    // read input image
    int w, h, pd;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
    double *out = malloc(w*h*pd*sizeof*out);

    // channel equalization
    channel_equalization(in, out, w, h, pd);

    // write output image
    iio_write_image_double_split(filename_out, out, w, h, pd);

    // free memory
    free(in);
    free(out);

    return EXIT_SUCCESS;
}
