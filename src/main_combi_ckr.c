// Irregularly sampled data fitting from images (mosaicked or not) using classical kernel regression

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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "xmtime.h"

#include "iio.h"
#include "combi_ckr_core.h"

int main(int c, char *v[])
{
    if (c != 10) {
            fprintf(stderr,"usage:\n\t%s path_image path_homo number_ini number_end zoom out order sigma raw\n", *v);
            //                         0 1          2         3          4          5    6   7     8     9
            fprintf(stderr,"example:\n\t%s im%%03d.tiff homo%%03d.tiff 1 30 1 output.tiff 2 0.7 0\n", *v);
            return EXIT_FAILURE;
    }

    // initialize time
    unsigned long t1 = xmtime();

    char *path_image = c > 1 ? v[1] : "-";
    char *path_homo = c > 2 ? v[2] : "-";
    int  ini = atoi(v[3]);
    int  end = atoi(v[4]);
    double zoom = atof(v[5]);
    char *filename_out= c > 6 ? v[6] : "-";
    int order = atoi(v[7]);
    float sigma2 = atof(v[8]);
    int raw = atoi(v[9]);
    
    if( sigma2 > 0 )
        sigma2 *= sigma2; // sigma*sigma
    else { // if sigma <= 0 a value of sigma is computed in order to have enough samples per pixel
        int Npoints = 25;
        int factorraw = 3*raw+(1-raw); // 3 if raw = 1 and 1 if raw = 0
        sigma2 = zoom*zoom*Npoints*factorraw/(3.14159*(end-ini+1)*4);
        sigma2 = (sigma2 < 0.04) ? 0.04 : sigma2; // sigma > 0.2
    }
    
    /* pre-computation and tabulation of the applicability function */
    double *applic= malloc((NPRECOMPUTE+1)*sizeof*applic);
    double sigma3 = R2MAX*1.0/(2*NPRECOMPUTE);
    for( int i=0; i<=NPRECOMPUTE; i++)
        applic[i] = exp(-i*sigma3); //exp(-x/(2*sigma2))

    int Nreg;
    switch( order ) {
        case 0:
        Nreg = 2+1+1;
        break;
        case 1:
        Nreg = 9+1+1;
        break;
        case 2:
        Nreg = 21+1+1;
        break;
        default:
        order = 2;
        Nreg  = 21+1+1;
        break;
    }
    
    char filename_image[500], filename_homo[500], filename_image2[500], filename_homo2[500];
    int w, h, pd, wout, hout;
    
    sprintf(filename_image, path_image, ini);
    double *image = iio_read_image_double_split(filename_image, &w, &h, &pd);
    free(image);
    if( raw ==1 && pd !=1 ) {
        printf("Raw input should have one channel \n");
        return EXIT_FAILURE;
    }
    int pd3 = (raw == 1) ? 3 : pd; // if raw input then the output image has 3 channels
    wout = w*zoom;
    hout = h*zoom;
    
    /* construct the array with the data used for the combination */
    double *buffer = malloc(Nreg*wout*hout*pd3*sizeof*buffer);
    double *buffer2 = malloc(Nreg*wout*hout*pd3*sizeof*buffer2);
    for(int i=0;i<Nreg*wout*hout*pd3; i++) {
        buffer[i] = 0;
        buffer2[i] = 0;
    }
    
    int number_total = end-ini+1;
    int number_set   = number_total/2;
    
    /* loop over the images */
    for (int ind=0; ind<number_set; ind++)
    {
        //printf("Computation for %i images\n",2*(ind+1));
        /* set 1 */
        sprintf(filename_image, path_image, ini+2*ind);
        sprintf(filename_homo, path_homo, ini+2*ind);
        
        /* set 2 */
        sprintf(filename_image2, path_image, ini+2*ind+1);
        sprintf(filename_homo2, path_homo, ini+2*ind+1);
        
        fill_buffer_routine_parallel(filename_image, filename_homo,filename_image2, filename_homo2,
                                     w, h, pd, buffer, buffer2, zoom, order, Nreg, sigma2, raw, applic);
    }

    if(number_total%2 == 1) { // odd number of input images
        sprintf(filename_image, path_image, end);
        sprintf(filename_homo, path_homo, end);
        fill_buffer_routine(filename_image, filename_homo, w, h, pd, buffer, zoom, order, Nreg, sigma2, raw, applic);
    }
    
    /* add the two buffers */
    for(int i=0;i<Nreg*wout*hout*pd3; i++)
        buffer[i] += buffer2[i];

    // final time and print time
    unsigned long t2 = xmtime();
    printf("Computation of matrice coefficients in %.3f seconds \n", (float) (t2-t1)/1000);
    
    // image computation
    double *out = malloc(wout*hout*pd3*sizeof*out);
    compute_image(out, buffer, wout, hout, pd3, order, Nreg);
    iio_write_image_double_split(filename_out, out, wout, hout, pd3);

    unsigned long t3 = xmtime();
    printf("Image computed in %.3f seconds \n", (float) (t3-t2)/1000);
    
    // free memory
    free(buffer);
    free(buffer2);
    free(out);
    free(applic);
    
    return EXIT_SUCCESS;
}
