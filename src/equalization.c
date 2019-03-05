// Program for equalizing the dynamics of two images (possibly mosaicked)
// Images should have the same number of channels
// The variants are:
// histogram equalization (histo) with a simple sort of the values (for images of the same size)
// affine equalization (affine) where the mean and the variance are equalized
// additive mean equalization (meanp)
// multiplicative mean equalization (meanx)

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

#include <string.h>

#include "equalization_core.c"
#include "iio.h"

int main(int c, char *v[])
{
  if (c < 5) {
    fprintf(stderr,"usage:\n\t%s ref modified out type [raw]\n", *v);
    //                         0 1   2        3   4    5
    fprintf(stderr,"type possibilities: histo, affine, meanp, meanx\n");
    fprintf(stderr,"Set raw=1 to handle RAW images\n");
    return EXIT_FAILURE;
  }

  char *filename_ref = c > 1 ? v[1] : "-";
  char *filename_modified = c > 2 ? v[2] : "-";
  char *filename_out = c > 3 ? v[3] : "-";
  char *equalization_type = c > 4 ? v[4] : "-";
  int raw = c > 5 ? atoi(v[5]) : 0;

  // read input images
  int w, h, pd;
  double *ref = iio_read_image_double_vec(filename_ref, &w, &h, &pd);

  int w2, h2 , pd2;
  double *modified = iio_read_image_double_vec(filename_modified, &w2, &h2, &pd2);

  // sanity checks
  if ( pd != pd2 ) {
    fprintf(stderr,"Images should have the same number of channel\n");
    return EXIT_FAILURE;
  }
  if( raw == 1) {
    if ( pd != 1 ) {
      fprintf(stderr,"Raw images must only have one channel\n");
      return EXIT_FAILURE;
    }
  }

  // equalization
  double *out = malloc(w2*h2*pd2*sizeof*out);
  if ( 0 == strcmp(equalization_type, "histo") ) {
    if( w!=w2 || h!=h2 || pd!=pd2 ) {
      fprintf(stderr,"Size of images should be the same\n");
      return EXIT_FAILURE;
    }
    equalization_histo(ref, modified, out, w, h, pd);
  }
  else if ( 0 == strcmp(equalization_type, "affine") ) {
    equalization_affine(ref, modified, out, w, h, w2, h2, pd, raw);
  }
  else if ( 0 == strcmp(equalization_type, "meanp") ) {
    equalization_meanp(ref, modified, out, w, h, w2, h2, pd, raw);
  }
  else if ( 0 == strcmp(equalization_type, "meanx") ) {
    equalization_meanx(ref, modified, out, w, h, w2, h2, pd, raw);
  }
  else {
    fprintf(stderr,"Unknown equalization type\n");
    return EXIT_FAILURE;
  }

  // write output image
  iio_write_image_double_vec(filename_out, out, w2, h2, pd2);

  // free memory
  free(ref);
  free(modified);
  free(out);

  return EXIT_SUCCESS;
}
