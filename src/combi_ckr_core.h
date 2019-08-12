#ifndef COMBI_CKR_ROUTINE_H
#define COMBI_CKR_ROUTINE_H

#define NPRECOMPUTE 1000000
#define R2MAX 16 // compute the weight for r2 in [0 rmax]

void compute_image(double *out, const double *buffer, int w, int h, int pd, int order, int Nreg);
void fill_buffer_routine(const char *filename_image, char *filename_homo, int w, int h, int pd,
                         double *buffer, float zoom, int order, int Nreg, double sigma2,
                         int raw, const double *applic);
void fill_buffer_routine_parallel(const char *filename_image, char *filename_homo,
                         const char *filename_image2, char *filename_homo2,
                         int w, int h, int pd,
                         double *buffer, double *buffer2, float zoom, int order, int Nreg, double sigma2,
                         int raw, const double *applic);


#endif
