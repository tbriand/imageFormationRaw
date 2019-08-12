#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"
#include <assert.h>

// We use the Bayer RGGB pattern
// The color channel in RAW images corresponds to i%2 + j%2
// R --> 0
// G --> 1
// B --> 2
static void rgb2raw(double *in, double *out, int w, int h)
{
    int canal;

    for(int j = 0; j < h; j++) {
        for(int i = 0; i < w; i++) {
            canal = i%2 + j%2;
            out[i + j*w] = in[i + j*w + canal*w*h];
        }
    }
}

int main(int c, char *v[])
{
    if (c != 3) {
            fprintf(stderr,"usage:\n\t%s in out\n", *v);
            //                         0 1  2    
            return EXIT_FAILURE;
    }
    char *filename_in = c > 1 ? v[1] : "-";
    char *filename_out = c > 2 ? v[2] : "-";
    
    // read input
    int w, h, pd;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &pd);
    if (pd != 3) {
        printf("Input image should be RGB \n");
        return EXIT_FAILURE;
    }
    
    double *out = malloc(w*h*sizeof*out);
    
    // compute the RAW version of the RGB input
    rgb2raw(in, out, w, h);
    
    // write the output
    iio_write_image_double_split(filename_out, out, w, h, 1);
    
    // free memory
    free(in);
    free(out);

    return EXIT_SUCCESS;
}
