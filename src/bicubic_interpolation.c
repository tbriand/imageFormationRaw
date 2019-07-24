#ifndef _BICUBIC_C
#define _BICUBIC_C

#include "getpixel.c"
#include <stdio.h>
#include <stdlib.h>

#ifndef BOUNDARY_DEFINITION
#define BOUNDARY_DEFINITION
typedef enum
{
    BOUNDARY_CONSTANT = 0,
    BOUNDARY_HSYMMETRIC = 1,
    BOUNDARY_WSYMMETRIC = 2,
    BOUNDARY_PERIODIC = 3
} BoundaryExt;
#endif

static double cubic_interpolation(double v[4], float x)
{
    return v[1] + 0.5 * x*(v[2] - v[0]
                    + x*(2.0*v[0] - 5.0*v[1] + 4.0*v[2] - v[3]
                    + x*(3.0*(v[1] - v[2]) + v[3] - v[0])));
}

static double bicubic_interpolation_cell(double p[4][4], float x, float y)
{
    double v[4];
    v[0] = cubic_interpolation(p[0], y);
    v[1] = cubic_interpolation(p[1], y);
    v[2] = cubic_interpolation(p[2], y);
    v[3] = cubic_interpolation(p[3], y);
    return cubic_interpolation(v, x);
}

void interpolate_bicubic(double *out, double *in, int w, int h, BoundaryExt bc, double *xpos, double *ypos, int numPixels) {

    // boundary handling
    getsample_operator_double p;
    if ( bc == BOUNDARY_PERIODIC )
        p = getsample_per_double;
    if ( bc == BOUNDARY_CONSTANT )
        p = getsample_constant_double;
    if ( bc == BOUNDARY_HSYMMETRIC )
        p = getsample_2_double;
    if ( bc == BOUNDARY_WSYMMETRIC )
        p = getsample_3_double;

    double x, y;
    for (int k = 0; k < numPixels; k++) {
        x = xpos[k] - 1;
        y = ypos[k] - 1;

        int ix = floor(x);
        int iy = floor(y);
        double c[4][4];
        for (int j = 0; j < 4; j++)
            for (int i = 0; i < 4; i++)
                c[i][j] = p(in, w, h, 1, ix + i, iy + j, 0);
        out[k] = bicubic_interpolation_cell(c, x - ix, y - iy);
    }
}

#endif//_BICUBIC_C
