// Concatenation of 4 images (of the same size) into a single one
// The common width is assumed to be even (otherwise a crop is made)

#include <stdio.h>
#include <stdlib.h>
#include "iio.h"
#include <math.h>

int main(int c, char *v[])
{
  if (c != 6) {
    fprintf(stderr, "usage:\n\t%s im1 im2 im3 im4 out\n", *v);
    //                          0 1   2   3   4   5
    return EXIT_FAILURE;
  }

  char *filename_in1 = c > 1 ? v[1] : "-";
  char *filename_in2 = c > 2 ? v[2] : "-";
  char *filename_in3 = c > 3 ? v[3] : "-";
  char *filename_in4 = c > 4 ? v[4] : "-";
  char *filename_out = c > 5 ? v[5] : "-";

  int w, h, pd, w2, h2, pd2;
  double *image_in1 = iio_read_image_double_vec(filename_in1, &w, &h, &pd);
  double *image_in2 = iio_read_image_double_vec(filename_in2, &w2, &h2, &pd2);
  if( w2 != w || h2 != h || pd2 != pd ) {
      fprintf(stderr,"Incompatible image size\n");
      return EXIT_FAILURE;
  }
  double *image_in3 = iio_read_image_double_vec(filename_in3, &w2, &h2, &pd2);
  if( w2 != w || h2 != h || pd2 != pd ) {
      fprintf(stderr,"Incompatible image size\n");
      return EXIT_FAILURE;
  }
  double *image_in4 = iio_read_image_double_vec(filename_in4, &w2, &h2, &pd2);
  if( w2 != w || h2 != h || pd2 != pd ) {
      fprintf(stderr,"Incompatible image size\n");
      return EXIT_FAILURE;
  }
  
  // to preserve the bayer aspect we need the width to be even
  int wout = 2*(w/2), hout = h;
  double *out = malloc(4*wout*hout*pd*sizeof*out);
  for (int j = 0; j<hout; j++) {
      for (int i = 0; i<wout; i++) {
          for( int l = 0; l<pd; l++) {
              out[ l + pd*(i+4*wout*j)] = image_in1[l + pd*(i+w*j)];
              out[ l + pd*(i+wout+4*wout*j)] = image_in2[l + pd*(i+w*j)];
              out[ l + pd*(i+2*wout+4*wout*j)] = image_in3[l + pd*(i+w*j)];
              out[ l + pd*(i+3*wout+4*wout*j)] = image_in4[l + pd*(i+w*j)];
          }
      }
  }
  iio_write_image_double_vec(filename_out, out, 4*wout, hout, pd);
  
  // free memory
  free(image_in1);
  free(image_in2);
  free(image_in3);
  free(image_in4);
  free(out);
  
  return EXIT_SUCCESS;
}
