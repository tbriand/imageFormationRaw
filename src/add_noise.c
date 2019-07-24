#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"
#include "xmtime.h"
#include "random.c"

int main(int c, char *v[])
{
    if (c != 4) {
            fprintf(stderr,"usage:\n\t%s sigma in out\n", *v);
            //                         0 1     2  3
            return EXIT_FAILURE;
    }
    float sigma = atof(v[1]);
    char *filename_in = c > 2 ? v[2] : "-";
    char *filename_out= c > 3 ? v[3] : "-";

    int w, h, pd;
    double *in = iio_read_image_double_vec(filename_in, &w, &h, &pd);

    // initialize seed
    xsrand(xmtime());

    // add noise if sigma > 0
    if( sigma > 0)
        for(int i = 0; i < w*h*pd; i++ )
            in[i] += sigma*random_normal();

    // save output
    iio_write_image_double_vec(filename_out, in, w, h, pd);

    // free memory
    free(in);

    return EXIT_SUCCESS;
}
