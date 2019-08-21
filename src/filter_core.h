#ifndef FILTER_CORE_H
#define FILTER_CORE_H

void apply_asymptotic_nc_filter_global(double *in, double *out, int w, int h, int pd, int order, double sigma, int inverse);
void apply_bilinear_luminance(double *in, double *out, int w, int h);
void apply_hvs_lowpass(double *double_in, double *double_out, int w, int h);
void apply_perfect_lowpass(double *in, double *out, int w, int h, int pd);
void apply_integration_kernel(double *in, double *out, int w, int h, int pd, double zoom);

#endif