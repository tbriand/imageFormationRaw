// utility functions regarding the application of filters

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

#ifndef FILTER_CORE
#define FILTER_CORE

#include "fft_core.c"

void multiply_filter(double *outhat, double *inhat, double *filter, int w, int h, int pd)
{
    double (*outhat2)[w][2*pd] = (void*) outhat;
    double (*inhat2)[w][2*pd] = (void*) inhat;
    double (*filter2)[w] = (void*) filter;
    
    double tmp;
    
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++)
        {
            tmp = filter2[j][i];
            for (int l = 0; l < 2*pd; l++)
                outhat2[j][i][l]=inhat2[j][i][l]*tmp;
        }
}

void compute_asymptotic_nc_filter(double *filter, int order, double sigma, int w, int h)
{
    double (*filter2)[w] = (void*) filter;

    int ind_w = (w+1)/2;
    int ind_h = (h+1)/2;
    double fw = 2*M_PI/((double) w);
    double fh = 2*M_PI/((double) h);
    double sigma2 = sigma*sigma*0.5;
    double x, y;
    
    if ( order <= 1 ) {
        for (int j = 0; j < h; j++)
        {
        y = (j<ind_h) ? j : j-h;
        y *= fh;
            for (int i = 0; i < w; i++)
            {
                x = (i<ind_w) ? i: i-w;
                x *= fw;
                x = x*x+y*y;
                filter2[j][i] = exp(-sigma2*x);
            }
        }
    }
    else { //order = 2
        for (int j = 0; j < h; j++)
        {
        y = (j<ind_h) ? j : j-h;
        y *= fh;
            for (int i = 0; i < w; i++)
            {
                x = (i<ind_w) ? i: i-w;
                x *= fw;
                x = x*x+y*y;
                filter2[j][i] = exp(-sigma2*x)*(1+sigma2*x);
            }
        }
    }
}

void apply_asymptotic_nc_filter(double *in, double *out, int w, int h, int pd, int order, double sigma, int inverse)
{
  double *inhat = xmalloc(w*h*2*pd*sizeof*inhat);
  double *outhat = xmalloc(w*h*2*pd*sizeof*outhat);
  double *filter = xmalloc(w*h*sizeof*filter);
  
  //compute input DFT
  fft_direct(inhat, in, w, h, pd);

  //compute filter
  compute_asymptotic_nc_filter(filter, order, sigma, w, h);
  
  //inverse if needed
  if( inverse == 1 ) {
      for (int i=0; i<w*h; i++)
          filter[i] = 1.0/filter[i];
  }
      
  //multiply filter
  multiply_filter(outhat, inhat, filter, w, h, pd);

  // FFT inverse
  fft_inverse(out, outhat, w, h, 2*pd);
  
  free(inhat);
  free(outhat);
}

void apply_asymptotic_nc_filter_global(double *in, double *out, int w, int h, int pd, int order, double sigma, int inverse)
{
  double *in_sym = xmalloc(4*w*h*pd*sizeof*in_sym);  
  symmetrization(in, in_sym, w, h, pd);
  
  apply_asymptotic_nc_filter(in_sym, in_sym, 2*w, 2*h, pd, order, sigma, inverse);
  
  double (*out2)[w][pd] = (void*) out;
  double (*in_sym2)[2*w][pd] = (void*) in_sym;
  for(int l=0;l<pd;l++)
    for(int i=0;i<w;i++)
      for(int j=0; j<h;j++)
        out2[j][i][l] = in_sym2[j][i][l];

  free(in_sym);
}

void bilinear_demosaicing(double *in, double *out, int w, int h) {
    //filters for bilinear interp
    int RB[3][3] = {
        {1,2,1},
        {2,4,2},
        {1,2,1}
    };
    int G[3][3] = {
        {0,1,0},
        {1,4,1},
        {0,1,0}
    };

    //CFA image with 0 when unknown
    double *CFA = malloc(3*w*h*sizeof*CFA);
    FORI(3*w*h)
        CFA[i] = 0;
    int l;
    FORI(w) {
        FORJ(h) {
            l = i%2+j%2;
            CFA[(j*w+i)*3+l] = in[j*w+i];
        }
    }

    //bilinear interp
    double tmpR, tmpG, tmpB;
    FORI(w) {
        FORJ(h) {
            tmpR = 0;
            tmpG = 0;
            tmpB = 0;
            for (int x=-1; x<=1 ; x++) {
                for (int y=-1; y<=1 ; y++) {
                    int indx=j+x;
                    int indy=i+y;
                    if( indx>=0 && indx<h && indy>=0 && indy<w ) {
                        tmpR += CFA[(indx*w+indy)*3]*RB[y+1][x+1];
                        tmpG += CFA[(indx*w+indy)*3+1]*G[y+1][x+1];
                        tmpB += CFA[(indx*w+indy)*3+2]*RB[y+1][x+1];
                    }
                }
            }
            out[(j*w+i)*3] = tmpR*0.25;
            out[(j*w+i)*3+1] = tmpG*0.25;
            out[(j*w+i)*3+2] = tmpB*0.25;
        }
    }
    
    //free memory
    free(CFA);
}

void apply_bilinear_luminance( double *in, double *out, int w, int h ) {
    
    double *rgb = malloc(3*w*h*sizeof*rgb);
    bilinear_demosaicing(in, rgb, w, h);
    
    FORI(w*h)
        out[i] = 0.3333333333333*(rgb[3*i] + rgb[3*i+1] + rgb[3*i+2]);
    
    free(rgb);
}

void apply_hvs_lowpass( double *double_in, double *double_out, int w, int h)
{
  double (*double_in2)[w] = (void*) double_in;
  double (*double_out2)[w] = (void*) double_out;

  int filter[11][11] = {
    {0,0,0,0,1,0,1,0,0,0,0},
    {0,0,0,-1,0,-2,0,-1,0,0,0},
    {0,0,1,1,2,1,2,1,1,0,0},
    {0,-1,1,-5,3,-9,3,-5,1,-1,0},
    {1,0,2,3,1,7,1,3,2,0,1},
    {0,-2,1,-9,7,104,7,-9,1,-2,0},
    {1,0,2,3,1,7,1,3,2,0,1},
    {0,-1,1,-5,3,-9,3,-5,1,-1,0},
    {0,0,1,1,2,1,2,1,1,0,0},
    {0,0,0,-1,0,-2,0,-1,0,0,0},
    {0,0,0,0,1,0,1,0,0,0,0}
  };

  double ifactor = 1.0/128;
  FORI(w){
    FORJ(h){
      double tmp=0;
      for (int x=-5; x<=5 ; x++)
      {
	for (int y=-5; y<=5 ; y++)
	{
	  int indx=j+x;
	  int indy=i+y;
	  if( indx>=0 && indx<h && indy>=0 && indy<w )
	    tmp+= double_in2[indx][indy]*filter[x+5][y+5];
	}
      }
      double_out2[j][i] = tmp*ifactor;
    }
  }
}

#endif