// Application of the VST from ponomarenko estimation

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

#include "ponomarenko_core.c"

int main(int c, char *v[])
{
    if (c != 8) {
            fprintf(stderr, "usage:\n\t%s inpat_file inpat_img input_ext outpat_img type ini end\n", *v);
                            //          0 1          2         3         4          5    6   7
            fprintf(stderr, "Ponomarenko noise estimations should be in files inpat_file_i.txt (for i=1:4)\n");
            return EXIT_FAILURE;
    }
    
    char *inpat_file = v[1];
    char *inpat_img = v[2];
    char *ext = v[3];
    char *outpat_img = v[4];
    int type = atoi(v[5]);
    int ini = atoi(v[6]);
    int end = atoi(v[7]);

    if ( type == 0 ) {
        double R[2], G[2], B[2];
        int index = -1;

        // regression on the 3 channels to compute the coefficient of the vst
        regression_3channels(inpat_file, index, R, G, B);
        
        if( R[0] < 0)
            R[0] = 0;
        if( G[0] < 0)
            G[0] = 0;
        if( B[0] < 0)
            B[0] = 0;
        
        printf("R[0] = %lf \t R[1] = %lf\n",R[0],R[1]);
        printf("G[0] = %lf \t G[1] = %lf\n",G[0],G[1]);
        printf("B[0] = %lf \t B[1] = %lf\n",B[0],B[1]);
        
        // apply the vst on the 3 channels of the raw images
        for (index = ini; index <= end; index++) 
            vst_3channels(inpat_img, ext, outpat_img, index, R, G, B);
    }
    else {
        vst_trapeze_multiple(inpat_file, inpat_img, ext, outpat_img, ini, end);
    }

    return EXIT_SUCCESS;
}
