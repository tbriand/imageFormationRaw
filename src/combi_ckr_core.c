// Utility functions for the irregularly sampled data fitting using classical kernel regression

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "iio.h"
#include "homography_core.h"
#include "fft_core.h"
#include "filter_core.h"
#include "combi_ckr_core.h"

static double compute_pixel_value(const double *C, int order)
{
  double p=0;
  int test = 1; // test = 1 means no error while test = 0 means error --> decrease the order
  switch (order)
  {
    case 0: ;
      if( C[0] > 0 )
        p = C[1]/C[0];
      else
          if( C[2] > 0)
              p = C[3]/C[2];
      break;
    case 1: ;
      if( C[9] > 6 ) { // need 3 points at least but if less it might give bad results
        double det, m1, m2, m3;

        m1 = C[5]*C[3]-C[4]*C[4];
        m2 = -(C[5]*C[1]-C[2]*C[4]);
        m3 = C[4]*C[1]-C[2]*C[3];
        det = C[0]*m1 - C[1]*m2 + C[2]*m3;
        
        if( det > 10e-15 ) {
          p = m1*C[6]+m2*C[7]+m3*C[8];
          p /= det;
        }
        else
          test = 0;
      }
      else
        test = 0;

      if ( test == 0 || p != p ) {
        /* if the pixel value can not be computed at order 1 switch to order 0 */
        double C2[4];
        C2[0] = C[0];
        C2[1] = C[6];
        C2[2] = C[9];
        C2[3] = C[10];
        p = compute_pixel_value(C2, 0);
      }
      break;
    case 2: ;
      if ( C[21] > 12 ) { // need 6 points at least but if less it might give bad results
        // construct the matrix A
        double A[36];

        A[35] = C[0];
        A[34] = C[1];
        A[33] = C[2];
        A[32] = C[3];
        A[31] = C[4];
        A[30] = C[5];

        A[29] = C[1];
        A[28] = C[3];
        A[27] = C[4];
        A[26] = C[6];
        A[25]= C[7];
        A[24]= C[8];

        A[23]= C[2];
        A[22]= C[4];
        A[21]= C[5];
        A[20]= C[7];
        A[19]= C[8];
        A[18]= C[9];

        A[17]= C[3];
        A[16]= C[6];
        A[15]= C[7];
        A[14]= C[10];
        A[13]= C[11];
        A[12]= C[12];

        A[11]= C[4];
        A[10]= C[7];
        A[9]= C[8];
        A[8]= C[11];
        A[7]= C[12];
        A[6]= C[13];

        A[5]= C[5];
        A[4]= C[8];
        A[3]= C[9];
        A[2]= C[12];
        A[1]= C[13];
        A[0]= C[14];

        // compute B
        double B[6];
        for (int i=0; i<6;i++) B[5-i]=C[15+i];

        // cholesky decomposition A=L*L' where L is lower triangular
        double L[36]={0};
        double factor = sqrt(A[0]);
        if( factor > 0 )
	  L[0] = factor;
	else
	  test = 0;

	for (int i=1; i<6; i++)
	  L[6*i]= A[6*i]/factor;
	for (int j=1; j<6; j++) {
	  double sum=0;
	  for (int k=0;k<j;k++) sum+= L[j*6+k]*L[j*6+k];
	  factor = sqrt(A[6*j+j] - sum);
	  if( factor > 0) {
	    double factor2 =  1/factor;
	    L[6*j+j] = factor;
	    for (int i=j+1; i<6; i++) {
	      double sum2=0;
	      for (int k=0;k<j;k++) sum2+=L[6*i+k]*L[6*j+k];
	      L[6*i+j]= (A[6*i+j] - sum2)*factor2;
	    }
	  }
	  else
	    test = 0;
	}

	// the system is B=L*L'*x
	// we first compute y such that B=L*y
	double y[6]={0};
	if( L[0] > 0 )
	    y[0]=B[0]/L[0];
	else
	  test = 0;

	for (int i=1;i<6;i++)
	{
	  double sum = 0;
	  for (int j=0;j<i;j++) sum += L[6*i+j]*y[j];
	  if( L[6*i+i]> 0 )
	    y[i]=(B[i]-sum)/L[6*i+i];
	  else
	    test = 0;
	}
	
	// then y=L'*x with L' upper triangular so the coefficient x[5] is computed with a simple division
	if( L[35] > 0)
	  p = y[5]/L[35];
	else
	  test = 0;
      }
      else
        test = 0;

      if ( test == 0 || p != p ) {
          /* if the pixel value can not be computed at order 2 switch to order 0*/
            double C0[4];
            C0[0] = C[0];
            C0[1] = C[15];
            C0[2] = C[21];
            C0[3] = C[22];
            p = compute_pixel_value(C0, 0);
      }
      break;
  }
  return(p);
}

void compute_image(double *out, const double *buffer, int w, int h, int pd, int order, int Nreg)
{
    double C[Nreg];
    for(int i=0; i<w*h*pd; i++) {
        for(int j=0; j<Nreg; j++)
            C[j] = buffer[Nreg*i+j];

        out[i] = compute_pixel_value(C, order);
    }
}

static void fill_buffer(double *buffer, const double *image, double *homo, float zoom, int w, int h, int pd,
                        int order, int Nreg, double sigma2, int raw, const double *applic)
{
    int wout = w*zoom;
    int hout = h*zoom;

    // values for the indices
    double radius2 = R2MAX*sigma2;
    double step = NPRECOMPUTE*1.0/radius2;
    double indmax = sqrt(radius2);

    int pd3 = (raw == 1) ? 3 : pd;
    //double (*buffer2)[wout][pd3][Nreg] = (void *) buffer;

    // invert the homography
    double ihomo[9];
    invert_homography(ihomo,homo);

    int l, i, j, pd2, m, n, ind;
    double Im, Im2;
    double x, y, xc, yc, wc, r2, xc2, yc2, xc3, yc3;
    double p[2], q[2];
  
    /* loop over the pixel of the image */
    for (j=0; j<h ; j++) {
        p[1] = j;
        for (i=0; i<w; i++) {
            // compute the value of h^{-1}(i,j)
            p[0] = i;
            apply_homography(q, p, ihomo);
            x = zoom*q[0];
            y = zoom*q[1];
            
            // search the pixel for which (x,y) is in their neighborhood (circle of radius radius)
            for (m = x - indmax - 1; m <= x + indmax + 1; m++) // only look for relevant indices
                if(m>=0 && m<wout) // in the image
                    for (n = y - indmax - 1; n <= y + indmax+1; n++) // only look for relevant indices
                        if(n>=0 && n<hout) { // in the image
                            xc = x - m;
                            yc = y - n;
                            xc2 = xc*xc;
                            yc2 = yc*yc;
                            r2 = xc2+yc2;
                            if( r2<radius2 ) {
                                wc = applic[(int) round(r2*step)]; //r2 is almost i*R2MAX/NPRECOMPUTE = i/step
                                for(l = 0; l<pd; l++) {
                                    pd2 = (raw == 1) ? i%2 + j%2 : l; // Bayer filter RGGB
                                    ind = (m + n*wout + pd2*wout*hout)*Nreg; // index in buffer
                                    Im = image[i+j*w+l*w*h];
                                    Im2 = wc*Im;

                                    switch ( order ) {
                                        case 0:
                                            buffer[ind + 0] += wc;
                                            buffer[ind + 1] += Im2;
                                            buffer[ind + 2] += 1;
                                            buffer[ind + 3] += Im;
                                            break;
                                        case 1:
                                            buffer[ind + 0] += wc;
                                            buffer[ind + 1] += wc*xc;
                                            buffer[ind + 2] += wc*yc;
                                            buffer[ind + 3] += wc*xc2;
                                            buffer[ind + 4] += wc*xc*yc;
                                            buffer[ind + 5] += wc*yc2;
                                            buffer[ind + 6] += Im2;
                                            buffer[ind + 7] += Im2*xc;
                                            buffer[ind + 8] += Im2*yc;
                                            buffer[ind + 9] += 1;
                                            buffer[ind + 10] += Im;
                                            // C[0] <--> A11 <--> sum w_i
                                            // C[1] <--> A12 <--> sum w_i x_i
                                            // C[2] <--> A13 <--> sum w_i y_i
                                            // C[3] <--> A22 <--> sum w_i x_i^2
                                            // C[4] <--> A23 <--> sum w_i x_i y_i
                                            // C[5] <--> A33 <--> sum w_i y_i^2
                                            // C[6] <--> B1 <--> sum w_i I_i
                                            // C[7] <--> B2 <--> sum w_i I_i x_i
                                            // C[8] <--> B3 <--> sum w_i I_i y_i
                                            break;
                                        case 2:
                                            xc3=xc2*xc;
                                            yc3=yc2*yc;
                                            buffer[ind + 0] += wc;
                                            buffer[ind + 1] += wc*xc;
                                            buffer[ind + 2] += wc*yc;
                                            buffer[ind + 3] += wc*xc2;
                                            buffer[ind + 4] += wc*xc*yc;
                                            buffer[ind + 5] += wc*yc*yc;
                                            buffer[ind + 6] += wc*xc3;
                                            buffer[ind + 7] += wc*xc2*yc;
                                            buffer[ind + 8] += wc*xc*yc2;
                                            buffer[ind + 9] += wc*yc3;
                                            buffer[ind + 10] += wc*xc2*xc2;
                                            buffer[ind + 11] += wc*xc3*yc;
                                            buffer[ind + 12] += wc*xc2*yc2;
                                            buffer[ind + 13] += wc*xc*yc3;
                                            buffer[ind + 14] += wc*yc2*yc2;
                                            buffer[ind + 15] += Im2;
                                            buffer[ind + 16] += xc*Im2;
                                            buffer[ind + 17] += yc*Im2;
                                            buffer[ind + 18] += xc2*Im2;
                                            buffer[ind + 19] += xc*yc*Im2;
                                            buffer[ind + 20] += yc2*Im2;
                                            buffer[ind + 21] += 1;
                                            buffer[ind + 22] += Im;
                                            // C[0] <--> A11 <--> sum w_i
                                            // C[1] <--> A12 <--> sum w_i x_i
                                            // C[2] <--> A13 <--> sum w_i y_i
                                            // C[3] <--> A14 = A22 <--> sum w_i x_i^2
                                            // C[4] <--> A15 = A23 <--> sum w_i x_i y_i
                                            // C[5] <--> A16 = A33 <--> sum w_i y_i^2
                                            // C[6] <--> A24 <--> sum w_i x_i^3
                                            // C[7] <--> A25 = A34 <--> sum w_i x_i^2 y_i
                                            // C[8] <--> A26 = A35 <--> sum w_i x_i y_i^2
                                            // C[9] <--> A36 <--> sum w_i y_i^3
                                            // C[10]<--> A44 <--> sum w_i x_i^4
                                            // C[11]<--> A45 <--> sum w_i x_i^3 y_i
                                            // C[12]<--> A46 = A55 <--> sum w_i x_i^2 y_i^2
                                            // C[13]<--> A56 <--> sum w_i x_i y_i^3
                                            // C[14]<--> A66 <--> sum w_i y_i^4
                                            // C[15]<--> B1 <--> sum w_i I_i
                                            // C[16]<--> B2 <--> sum w_i I_i x_i
                                            // C[17]<--> B3 <--> sum w_i I_i y_i
                                            // C[18]<--> B4 <--> sum w_i I_i x_i^2
                                            // C[19]<--> B5 <--> sum w_i I_i x_i y_i
                                            // C[20]<--> B6 <--> sum w_i I_i y_i^2
                                            break;
                                        }
                                }
                            }
                        }
        }
    }
}

void fill_buffer_routine(const char *filename_image, char *filename_homo, int w, int h, int pd,
                         double *buffer, float zoom, int order, int Nreg, double sigma2,
                         int raw, const double *applic)
{
    /* read input image and verify the size */
    int w2, h2, pd2;
    double *image = iio_read_image_double_split(filename_image, &w2, &h2, &pd2);
    if( w!=w2 || h!=h2 || pd!=pd2 ) {
        printf("Size of images should be the same \n");
    }
    
    /* read input homography */
    double homo[9];
    read_homography_file(homo, filename_homo);
    
    /* fill buffer */
    fill_buffer(buffer, image, homo, zoom, w, h, pd, order, Nreg, sigma2, raw, applic);
    
    /* free memory */
    free(image);
}

void fill_buffer_routine_parallel(const char *filename_image, char *filename_homo,
                         const char *filename_image2, char *filename_homo2,
                         int w, int h, int pd,
                         double *buffer, double *buffer2, float zoom, int order, int Nreg, double sigma2,
                         int raw, const double *applic)
{
    /* read input image and verify the size */
    int w2, h2, pd2;
    double *image = iio_read_image_double_split(filename_image, &w2, &h2, &pd2);
    if( w!=w2 || h!=h2 || pd!=pd2 )
        printf("Size of images should be the same \n");
    double *image2 = iio_read_image_double_split(filename_image2, &w2, &h2, &pd2);
    if( w!=w2 || h!=h2 || pd!=pd2 )
        printf("Size of images should be the same \n");
    
    /* read input homography */
    double homo[9];
    read_homography_file(homo, filename_homo);
    
    double homo2[9];
    read_homography_file(homo2, filename_homo2);
    
    /* fill buffer */
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
            fill_buffer(buffer, image, homo, zoom, w, h, pd, order, Nreg, sigma2, raw, applic);
        }
        #pragma omp section
        {
            fill_buffer(buffer2, image2, homo2, zoom, w, h, pd, order, Nreg, sigma2, raw, applic);
        }
    }

    /* free memory */
    free(image);
    free(image2);
}

