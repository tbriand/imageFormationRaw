// Functions around homographies

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

#ifndef HOMOGRAPHY_CORE
#define HOMOGRAPHY_CORE

#include <assert.h>
#include <math.h>

#include "random.h"
#include "parsenumbers.c"
#include "xfopen.c"
#include "cmphomod.h"


static void convert_ica_homography(double *homo) {
    if( homo[0] == 2 ) {
        // (2, h2, h5)
            homo[5] = homo[2];
            homo[2] = homo[1];
            homo[0] = 1;
            homo[1] = 0;
            homo[3] = 0;
            homo[4] = 1;
            homo[6] = 0;
            homo[7] = 0;
            homo[8] = 1;
    }
    if( homo[0] == 3 ) {
        // (3, h2, h5, THETA)
            homo[5] = homo[2];
            homo[2] = homo[1];
            homo[0] = cos(homo[3]);
            homo[1] = -sin(homo[3]);
            homo[3] = -homo[1];
            homo[4] = homo[0];
            homo[6] = 0;
            homo[7] = 0;
            homo[8] = 1;
    }
    if( homo[0] == 4 ) {
        // (4, h2, h5, h0-1, -h1)
            homo[5] = homo[2];
            homo[2] = homo[1];
            homo[0] = homo[3] + 1;
            homo[1] = -homo[4];
            homo[3] = -homo[1];
            homo[4] = homo[0];
            homo[6] = 0;
            homo[7] = 0;
            homo[8] = 1;
    }
    if( homo[0] == 6 ) {
        // (6, h2, h5, h0-1, h1, h3, h4-1)
        double tx = homo[2];
        double ty = homo[1];
        homo[0] = homo[3] + 1;
        homo[1] = homo[4];
        homo[2] = ty;
        homo[3] = homo[5];
        homo[4] = homo[6] + 1;
        homo[5] = tx;
        homo[6] = 0;
        homo[7] = 0;
        homo[8] = 1;
//for (int i=0; i<9;i++) printf("homo[%i] = %lf\n",i, homo[i]);
    }
    if( homo[0] == 8 ) {
        // (8, h0-1, h1, h2, h3, h4-1, h5, h6, h7)
            homo[0] = homo[1] + 1;
            homo[1] = homo[2];
            homo[2] = homo[3];
            homo[3] = homo[4];
            homo[4] = homo[5] + 1;
            homo[5] = homo[6];
            homo[6] = homo[7];
            homo[7] = homo[8];
            homo[8] = 1;
    }
}

static void read_homography(double *homo, char *characters) {
    int maxparam = 9;
    int nparams = parse_doubles(homo, maxparam, characters);
    assert(nparams <= maxparam);
    if ( homo[0] + 1 == nparams )
        convert_ica_homography(homo);
}

static void read_homography_file(double *homo, char *filename) {
    FILE *f = xfopen(filename, "r");
    int maxparam = 9;
    int nparams;
    double *homo_tmp = read_ascii_doubles(f, &nparams);
    xfclose(f);
    assert(nparams <= maxparam);
    for (int i=0; i<nparams; i++)
        homo[i] = homo_tmp[i];
    free(homo_tmp);
    if ( homo[0] + 1 == nparams )
        convert_ica_homography(homo);
}

static void invert_homography(double iH[9], double H[9])
{
  double det = H[0] * (H[4]*H[8] - H[5] * H[7]);
  det -= H[1] * (H[3]*H[8] - H[5] * H[6]);
  det += H[2] * (H[3]*H[7] - H[4] * H[6]);

  double tmp = 1/det;

  iH[0] = tmp * (H[4] * H[8] - H[5] * H[7]);
  iH[3] = tmp * (H[5] * H[6] - H[3] * H[8]);
  iH[6] = tmp * (H[3] * H[7] - H[4] * H[6]);

  iH[1] = tmp * (H[2] * H[7] - H[1] * H[8]);
  iH[4] = tmp * (H[0] * H[8] - H[2] * H[6]);
  iH[7] = tmp * (H[1] * H[6] - H[0] * H[7]);

  iH[2] = tmp * (H[1] * H[5] - H[2] * H[4]);
  iH[5] = tmp * (H[2] * H[3] - H[0] * H[5]);
  iH[8] = tmp * (H[0] * H[4] - H[1] * H[3]);
}

static void apply_homo(double y[2], double x[2], double H[9])
{
  double z[3];
  z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
  z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
  z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
  y[0] = z[0]/z[2];
  y[1] = z[1]/z[2];
}

static void apply_homography(double y[2], double x[2], double H[9])
{
  double z[3];
  z[0] = H[0]*x[0] + H[1]*x[1] + H[2];
  z[1] = H[3]*x[0] + H[4]*x[1] + H[5];
  z[2] = H[6]*x[0] + H[7]*x[1] + H[8];
  y[0] = z[0]/z[2];
  y[1] = z[1]/z[2];
}

// compute homo3 = homo2 o homo1
static void compose_homography(double *homo1, double *homo2, double *homo3)
{
    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            homo3[3*i+j]=0;
            for(int k=0;k<3;k++) homo3[3*i+j] += homo2[3*i+k]*homo1[3*k+j];
        }
    }
}

// translate homography in order to compute the result of a crop
// h_in(x+tx,y+ty) = (tx,ty) + h_out(x,y)
// i.e h_out = t(-tx,-ty) o h_in o t(tx,ty)
static void translate_homography(double *h_in, double *h_out, double tx, double ty)
{
    h_out[8] = h_in[6]*tx + h_in[7]*ty + h_in[8];
    h_out[0] = h_in[0] - tx*h_in[6];
    h_out[1] = h_in[1] - tx*h_in[7];
    h_out[2] = (h_in[0]*tx + h_in[1]*ty + h_in[2]) - tx*h_out[8];
    h_out[3] = h_in[3] - ty*h_in[6];
    h_out[4] = h_in[4] - ty*h_in[7];
    h_out[5] = (h_in[3]*tx + h_in[4]*ty + h_in[5]) - ty*h_out[8];
    h_out[6] = h_in[6];
    h_out[7] = h_in[7];
}

// zoom homography in order to be compatible with a zoom of an image
// (zoomx, zoomy) h_in (x/zoomx,y/zoomy) = h_out(x,y)
static void zoom_homography(double *h_in, double *h_out, double zx, double zy)
{
    h_out[0] = h_in[0];
    h_out[1] = h_in[1]*zx*1.0/zy;
    h_out[2] = h_in[2]*zx;
    h_out[3] = h_in[3]*zy*1.0/zx;
    h_out[4] = h_in[4];
    h_out[5] = h_in[5]*zy;
    h_out[6] = h_in[6]/zx;
    h_out[7] = h_in[7]/zy;
    h_out[8] = h_in[8];
}

static void set_to_zero( double *in, int w, int h, int pd, double H[9]) {    
    // set to 0 the value of the array if extrapolated value
    double p[2], q[2];
    for( int j=0; j<h; j++) {
        p[1] = j;
        for( int i=0; i<w; i++) {
            p[0] = i;
            apply_homo(q, p, H);
            if( q[0] < 0 || q[0] > w-1 || q[1] < 0 || q[1] > h-1 ) {
                for( int l = 0; l<pd; l++)
                    in[l + pd*(i+j*w) ] = 0.0;
            }
        }
    }
}

static void set_to_zero_split( double *in, int w, int h, int pd, double H[9]) {
    // set to 0 the value of the array if extrapolated value
    double p[2], q[2];
    for( int j=0; j<h; j++) {
        p[1] = j;
        for( int i=0; i<w; i++) {
            p[0] = i;
            apply_homo(q, p, H);
            if( q[0] < 0 || q[0] > w-1 || q[1] < 0 || q[1] > h-1 ) {
                for( int l = 0; l<pd; l++)
                    in[i+j*w + l*w*h] = 0.0;
            }
        }
    }
}

// compute the field as norm(homo2^{-1}( homo1( x ) ) - x)
static void compute_field_inv(float *field, double *homo1, double *homo2, int w, int h)
{
    float (*field2)[w] = (void *) field;
    double x[2];
    double y[2];
    double z[2];
    double ihomo2[9];
  
    invert_homography(ihomo2, homo2);
    
    for( int i = 0; i < w; i++ ){
        x[0] = i;
        for( int j = 0; j < h; j++){
            // initialisation
            x[1] = j;
      
            // calcul de h2^-1 \circ h - Id
            apply_homo(y,x,homo1);
            apply_homo(z,y,ihomo2);
            z[0] -= x[0];
            z[1] -= x[1];
      
            // fill the output image
            field2[j][i] = hypot(z[0],z[1]);
        }
    }
}

// compute the field as norm(homo2(x) - homo1(x))
static void compute_field(float *field, double *homo1, double *homo2, int w, int h)
{
    float (*field2)[w] = (void *) field;
    double x[2];
    double z1[2];
    double z2[2];

    for( int i = 0; i < w; i++ ){
        x[0] = i;
        for( int j = 0; j < h; j++){
            // initialisation
            x[1] = j;
      
            // calcul de h2^-1 \circ h - Id
            apply_homo(z1,x,homo1);
            apply_homo(z2,x,homo2);
            
            // fill the output image
            field2[j][i] = hypot(z1[0]-z2[0],z1[1]-z2[1]);
        }
    }
}

static void create_random_translation(double H[9], int L)
{
        double tx = (2*L*random_uniform() - L);
        double ty = (2*L*random_uniform() - L);
        H[0] = 1;
        H[1] = 0;
        H[2] = tx;
        H[3] = 0;
        H[4] = 1;
        H[5] = ty;
        H[6] = 0;
        H[7] = 0;
        H[8] = 1;
}

static void create_random_homography(double H[9], int w, int h, int L)
{
    double corner[4][2] = {{0,0}, {w-1,0}, {0,h-1}, {w-1,h-1}};
    double corner2[4][2];
    double a;
    
    for(int j=0; j<2; j++)
        for(int i=0; i<4; i++) {
            a = random_uniform();
            corner2[i][j] = corner[i][j] + (2*L*a - L);
        }

    double R[3][3];
    homography_from_4corresp(
    corner[0], corner[1], corner[2], corner[3],
    corner2[0], corner2[1], corner2[2], corner2[3], R);

    for(int i=0; i<9; i++) 
        H[i] = R[i/3][i%3];
}

#endif
