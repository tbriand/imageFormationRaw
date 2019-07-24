#include "homography_core.c"
#include "interpolation_core.c"
#include "iio.h"

int main(int c, char *v[])
{
    if (c < 3) {
    fprintf(stderr,"usage:\n\t%s in base_out number [interp boundary L type zoom]  \n", *v);
    //                        0  1  2        3       4      5        6 7    8
    return EXIT_FAILURE;
    }

    char *filename_in = v[1];
    char *base_out = v[2] ;
    int n = atoi(v[3]);
    char *interp = c > 4 ? v[4] : "splineper";
    char *boundary = c > 5 ? v[5] : "hsym";
    int L = c > 6 ? atoi(v[6]) : 3;
    int type = c > 7 ? atoi(v[7]) : 8;
    float zoom = c > 8 ? atof(v[8]) : 1.0;

    // read data
    int w, h, pd;
    double *in = iio_read_image_double_vec(filename_in, &w, &h, &pd);

    // rand initialization
    xsrand(0);

    // create all the homographies
    // it is done in the beginning so that for a given seed we always have the same homographies
    double *homographies = malloc(9*n*sizeof(double));
    double H[9];
    char filename_homo[500];

    // randomly compute others
    for (int j=0; j<n;j++) {
        if ( j > 0 ) {
            if( type == 2 ) {
                create_random_translation(H, L);
            }
            else if ( type == 3 ) {
                create_random_translation(H, L);
                double theta = 5*M_PI/180*(2*random_uniform()-1); //rotation of at most 5°
                int w2 = (w+1)/2;
                int h2 = (h+1)/2;
                H[0] = cos(theta);
                H[1] = -sin(theta);
                H[3] = -H[1];
                H[4] = H[0];
                H[2] += w2*(1-H[0])-h2*H[1];
                H[5] += -w2*H[3]+h2*(1-H[4]);
            }
            else if ( type == 6 ) {
                create_random_homography(H, w, h, L);
                H[6] = 0;
                H[7] = 0;
            }
            else if ( type == 8 ) {
                create_random_homography(H, w, h, L);
            }
            else {
                fprintf(stderr,"Wrong type of transformations...\n");
                return EXIT_FAILURE;
            }
        }
        else { // fill first with identity
            H[0] = 1.0;
            H[1] = 0.0;
            H[2] = 0.0;
            H[3] = 0.0;
            H[4] = 1.0;
            H[5] = 0.0;
            H[6] = 0.0;
            H[7] = 0.0;
            H[8] = 1.0;
        }
        
        sprintf(filename_homo, "%s_%i.hom", base_out, j+1);
        FILE *f = xfopen(filename_homo, "w");
        for (int i = 0; i < 9; i++) {
            homographies[9*j+i] = H[i];
            fprintf(f,"%1.16lg%c", H[i], i==8?'\n':' ');
        }
        xfclose(f);
    }

    //Boudary condition 
    BoundaryExt boundaryExt = read_ext(boundary);  

    /* compute images */
    int wout = w/zoom;
    int hout = h/zoom;
    char filename_out[500];

    for (int j = 0; j < n; j++) {
        double *out = malloc(wout*hout*pd*sizeof*out);
        for(int i = 0; i < 9; i++)
            H[i] = homographies[9*j+i];
        
        resampling(out, in, w, h, pd, H, interp, boundaryExt, zoom);
        
        sprintf(filename_out, "%s_%i.tiff", base_out, j+1);
        iio_write_image_double_vec(filename_out, out, wout, hout, pd);
        free(out);
    }

    // free memory
    free(in);
    free(homographies);

    return EXIT_SUCCESS;
}
