#ifndef PONOMARENKO_CORE_H
#define PONOMARENKO_CORE_H

void regression_1channel(char *filename, double C[2]);
void regression_3channels(char *inpat, int ind, double R[2], double G[2], double B[2]);
void vst_1channel(char *filename_in, char *filename_out, const double C[2]);
void vst_3channels(char *inpat, char *ext, char *outpat, int ind, double R[2], double G[2], double B[2]);
void vst_trapeze(char *inpat, char *image_in, char *image_out);
void vst_trapeze_multiple(char *inpat, char *inpat_img, char *ext, char *outpat_img, int ini, int end);

#endif
