#ifndef BICUBIC_H
#define BICUBIC_H

void interpolate_bicubic(double *out, double *in, int w, int h, int pd, BoundaryExt bc, double *xpos, double *ypos, int numPixels);

#endif
