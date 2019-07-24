#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iio.h"
#include "homography_core.c"

int main(int c, char *v[])
{
    if (c < 6) {
            fprintf(stderr, "usage:\n\t%s w h \"homo 1\" \"homo 2\""
                            //         0  1 2 3          4 
                            " field opt\n", *v);
                            //5     6
            return EXIT_FAILURE;
    }
    
    int w = atoi(v[1]);
    int h = atoi(v[2]);
    int opt = c > 6 ? atoi(v[6]) : 0;
    
    float *field = malloc(w*h*sizeof*field);

    int maxparam = 9;
    double homo1[maxparam];
    double homo2[maxparam];
    
    read_homography(homo1, v[3]);
    read_homography(homo2, v[4]);
    
    if ( opt == 0 ) {
        compute_field_inv(field, homo1, homo2, w, h);
    }
    else {
        compute_field(field, homo1, homo2, w, h);
    }
    
    iio_write_image_float_vec(v[5], field, w, h, 1);

    free(field);

    return EXIT_SUCCESS;
}
