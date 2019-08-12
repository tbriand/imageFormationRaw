/* SPDX-License-Identifier: GPL-2.0+
 * 
 * Thibaud Briand <briand.thibaud@gmail.com>
 * 
 * Copyright (c) 2018-2019, Thibaud Briand
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "iio.h"
#include "xmtime.h"
#include "interpolation_core.h"
#include "homography_core.h"
#include "fft_core.h"

// Function to transform char of the form "v0 v1 ..." into an array
// of doubles t where t[i] = vi.
// Taken from parsenumbers.c in the imscript repository:
// https://github.com/mnhrdt/imscript
static int parse_doubles(double *t, int nmax, const char *s)
{
    int i = 0, w;
    while (i < nmax && 1 == sscanf(s, "%lg %n", t + i, &w)) {
            i += 1;
            s += w;
    }
    return i;
}

// display help usage
void print_help(char *name)
{
    printf("\n<Usage>: %s input output \"h11 h12 h13 h21 h22 h23 h31 h32 h33\" [interp bc inverse]\n\n", name);
}

// Main function for the geometric transformation of an image
// using the NFFT algorithm (Algorithm 2).
int main(int c, char *v[])
{
    // initialize FFTW
    init_fftw();

    if (c < 3) {
        print_help(v[0]);
        return EXIT_SUCCESS;
    }
    
    char *filename_in  = v[1];
    char *filename_out = v[2];
    char *input_params = v[3];
    char *interp = c > 4 ? v[4] : "p+s-spline11";
    char *boundary = c > 5 ? v[5] : "hsym";
    int inverse = c > 6 ? atoi(v[6]) : 0;
    
    // read image
    int w, h, pd;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
    
    // initialize time
    unsigned long t1 = xmtime();
    
    // Read transformation
    double H[9];
    int maxparam = 9;
    int nparams = parse_doubles(H, maxparam, input_params);
    if ( nparams != maxparam ) {
        fprintf(stderr,"Incorrect input homography\n");
        return EXIT_FAILURE;
    }
    
    // Inverse homography if asked
    if ( inverse ) {
        double iH[9];
        invert_homography(iH, H);
        memcpy(H, iH, 9*sizeof(double));
    }
    
    //Boudary condition 
    BoundaryExt boundaryExt = read_ext(boundary);
    
    // memory allocation
    double *out = malloc(w*h*pd*sizeof*out);
    
    // homographic transformation of the image
    interpolate_image_homography(out, in, w, h, pd, H, interp, boundaryExt, 1);

    // final time and print time
    unsigned long t2 = xmtime();
    printf("Interpolation made in %.3f seconds \n", (float) (t2-t1)/1000);

    // write output image
    iio_write_image_double_split(filename_out, out, w, h, pd);
    
    // free memory
    free(in);
    free(out);
    clean_fftw();
    
    return EXIT_SUCCESS;
}
