#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

int main(int c, char *v[])
{
    if (c != 4) {
            fprintf(stderr,"usage:\n\t%s im1   im2  out\n", *v);
            //                         0 1     2    3
            return EXIT_FAILURE;
    }

    char *filename_im1 = c > 1 ? v[1] : "-";
    char *filename_im2 = c > 2 ? v[2] : "-";
    char *filename_out = c > 3 ? v[3] : "-";

    int w, h, pd;
    double *im1 = iio_read_image_double_split(filename_im1, &w, &h, &pd);
    int w2, h2 , pd2;
    double *im2 = iio_read_image_double_split(filename_im2, &w2, &h2, &pd2);
    
    if( w!=w2 || h!=h2 || pd!=pd2 ) {
        printf("Size of images should be the same \n");
            return EXIT_FAILURE;
    }
    
    double *out = malloc(w*h*pd*sizeof*out);
    for( int i=0; i<pd*w*h; i++) 
            out[i]=im1[i]-im2[i];
    
    iio_write_image_double_split(filename_out, out, w, h, pd);

    free(im1);
    free(im2);
    free(out);
    
    return EXIT_SUCCESS;
} 
