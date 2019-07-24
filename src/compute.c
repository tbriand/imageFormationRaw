#include <stdio.h>
#include <string.h>
#include "iio.h"
#include "compute_core.c"

int main(int c, char *v[])
{
    // display usage
    if (c < 4) {
        fprintf(stderr,"usage:\n\t%s type centered im [canal [raw]]\n", *v);
        //                         0 1    2        3  4       5
        fprintf(stderr,"type: rmse, max, max_abs, min, min_abs, mean, mean_abs, std\n");
        fprintf(stderr,"centered: 1 --> centered, other value --> whole image,\n");
        fprintf(stderr,"canal: integer value between 1 and the number of input channels\n");
        fprintf(stderr,"raw : set to 1 for raw images. Images are supposed to have a RGGB pattern\n");
        return EXIT_FAILURE;
    }
    
    // read parameters
    char *type = c > 1 ? v[1] : "-";
    int centered = c > 2 ? atoi(v[2]) : 0;
    char *in = c > 3 ? v[3] : "-";
    int canal = c > 4 ? atoi(v[4]) : 0;
    int raw = c > 5 ? atoi(v[5]) : 0;
    
    // read image
    int w, h, pd;
    double *x = iio_read_image_double_vec(in, &w, &h, &pd);
    
    if ( raw == 1 && pd != 1 ) {
        fprintf(stderr,"Raw images should have one channel\n");
        return EXIT_FAILURE;
    }
    if ( (raw != 1 && canal > pd) || canal > 3) {
        fprintf(stderr,"Channel number is too high\n");
        return EXIT_FAILURE;
    }

    // compute indices
    int xi = (centered == 1) ? w/3 : 0;
    int xf = (centered == 1) ? 2*w/3 : w;
    int w2 = xf - xi;
    
    int yi = (centered == 1) ? h/3 : 0;
    int yf = (centered == 1) ? 2*h/3 : h;
    int h2 = yf - yi;
    
    int pd2 = ( canal > 0 ) ? 1 : pd;
    double *y = malloc(w2*h2*pd2*sizeof*y);

    // image of interest (canal + eventual centering)
    int N; 
    if ( canal > 0 ) { // if the channel is specified
        if ( raw == 1 ) {
            N = 0;
            int l;
            for(int j=yi; j<yf; j++) {
                for(int i=xi; i<xf; i++) {
                    l = i%2 + j%2 + 1;
                    if ( l == canal ) {
                        y[N] = x[i+w*j];
                        N++;
                    }
                }
            }
        }
        else {
            N = w2*h2;
            for(int j=0; j<h2; j++) {
                for(int i=0; i<w2; i++) {
                    y[i+w2*j] = x[canal-1+pd*(i+xi+w*(j+yi))];
                }
            }
        }
    }
    else { // if the channel is not specified
        N = w2*h2*pd2;
        for(int j=0; j<h2; j++) {
            for(int i=0; i<w2; i++) {
                for(int l=0; l<pd2; l++) {
                    y[l+pd2*(i+w2*j)] = x[l+pd*(i+xi+w*(j+yi))];
                }
            }
        }
    }
    
    // compute value
    double val;
    if (0 == strcmp(type, "rmse")) {
        val = rmse(y, N);
    }
    if (0 == strcmp(type, "std")) {
        val = mean(y, N);
        val = std(y, val, N);
    }
    else if (0 == strcmp(type, "max")) {
        val = max(y, N);
    }
    else if (0 == strcmp(type, "max_abs")) {
        for (int i=0; i<N; i++)
            y[i] = fabs(y[i]);
        val = max(y, N);
    }
        else if (0 == strcmp(type, "min")) {
        val = min(y, N);
    }
    else if (0 == strcmp(type, "min_abs")) {
        for (int i=0; i<N; i++)
            y[i] = fabs(y[i]);
        val = min(y, N);
    }
    else if (0 == strcmp(type, "mean")) {
        val = mean(y, N);
    }
    else if (0 == strcmp(type, "mean_abs")) {
        for (int i=0; i<N; i++)
            y[i] = fabs(y[i]);
        val = mean(y, N);
    }
    
    // print value
    printf("%1.14lg\n", val);
    
    free(x);
    free(y);
    return EXIT_SUCCESS;
}