#ifndef INTERPOLATION_CORE_H
#define INTERPOLATION_CORE_H

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

BoundaryExt read_ext(const char* boundary);
void interpolate_image_homography(double *out, double *in, int w, int h, int pd, double H[9], 
                                  char *interp, BoundaryExt boundaryExt, double zoom);

#endif