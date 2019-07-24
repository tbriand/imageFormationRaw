#ifndef _GETPIXEL_C2
#define _GETPIXEL_C2

typedef float (*getsample_operator)(float*,int,int,int,int,int,int);
typedef double (*getsample_operator_double)(double*,int,int,int,int,int,int);
//typedef void (*setsample_operator)(float*,int,int,int,int,int,int,float);

// extrapolate by 0
inline
static float getsample_0(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return 0;
	return x[(i+j*w)*pd + l];
}

// extrapolate by nearest value
inline
static float getsample_1(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (l < 0) l = 0;
	if (i >= w) i = w-1;
	if (j >= h) j = h-1;
	if (l >= pd) l = pd-1;
	return x[(i+j*w)*pd + l];
}

#ifdef NAN
// extrapolate by nan
inline
static float getsample_nan(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return NAN;
	return x[(i+j*w)*pd + l];
}
#endif//NAN

// force a segfault if extrapolation is required
inline
static float getsample_error(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return *(volatile int*)0;
	return x[(i+j*w)*pd + l];
}

// abort the program extrapolation is required
inline
static float getsample_abort(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		abort();
	return x[(i+j*w)*pd + l];
}

// exit the program extrapolation is required
inline
static float getsample_exit(float *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		exit(42);
	return x[(i+j*w)*pd + l];
}

static int good_modulus(int n, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(n, -p);

	int r;
	if (n >= 0)
		r = n % p;
	else {
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	//assert(r >= 0);
	//assert(r < p);
	return r;
}

// static int positive_reflex(int n, int p)
// {
// 	int r = good_modulus(n, 2*p);
// 	if (r == p)
// 		r -= 1;
// 	if (r > p)
// 		r = 2*p - r;
// 	if (n < 0 && p > 1) r += 1;
// 	//assert(r >= 0);
// 	//assert(r < p);
//         if( r == p) {
//             printf("input n=%i p=%i\n",n,p);
//             printf("modulus %i\n",good_modulus(n, 2*p));
//             printf("result r=%i\n",r);
//         }
// 	return r;
// }

inline static int positive_reflex(int i, int N) {
    while(1) {
        if(i < 0)
            i = -1-i;
        else if(i >= N)
            i = (2*N-1)-i;
        else
            return i;
    }
}

// // ok for n in [-2p+1,2p-1]
// static int positive_reflex2(int n, int p)
// {
// 	int r = good_modulus(n, 2*p);
// 	if (r >= p && r<2*p-1)
// 		r = 2*p - r - 2;
// 	if (r==2*p-1)
// 		r = 1;
// 	//assert(r >= 0);
// 	//assert(r < p);
// 	return r;
// }

inline static int positive_reflex2(int i, int N) {
    while(1) {
        if(i < 0)
            i = -i;
        else if(i >= N)
            i = (2*N-2)-i;
        else
            return i;
    }
}

/* float */
// extrapolate by reflection (hsym)
inline
static float getsample_2(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = positive_reflex(i, w);
	j = positive_reflex(j, h);
	return getsample_abort(x, w, h, pd, i, j, l);
}

// extrapolate by reflection (wsym)
inline
static float getsample_3(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = positive_reflex2(i, w);
	j = positive_reflex2(j, h);
	return getsample_abort(x, w, h, pd, i, j, l);
}

// extrapolate by periodicity
inline
static float getsample_per(float *x, int w, int h, int pd, int i, int j, int l)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	return getsample_abort(x, w, h, pd, i, j, l);
}

// extrapolate by constant (set by calling it with zero sizes)
inline static
float getsample_constant(float *x, int w, int h, int pd, int i, int j, int l)
{
	static float value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return value;
	return x[(i+j*w)*pd + l];
}


/* double */

// abort the program extrapolation is required
inline
static double getsample_abort_double(double *x, int w, int h, int pd, int i, int j, int l)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		abort();
	return x[(i+j*w)*pd + l];
}

// extrapolate by reflection (hsym)
inline
static double getsample_2_double(double *x, int w, int h, int pd, int i, int j, int l)
{
	i = positive_reflex(i, w);
	j = positive_reflex(j, h);
	return getsample_abort_double(x, w, h, pd, i, j, l);
}

// extrapolate by reflection (wsym)
inline
static double getsample_3_double(double *x, int w, int h, int pd, int i, int j, int l)
{
	i = positive_reflex2(i, w);
	j = positive_reflex2(j, h);
	return getsample_abort_double(x, w, h, pd, i, j, l);
}

// extrapolate by periodicity
inline
static double getsample_per_double(double *x, int w, int h, int pd, int i, int j, int l)
{
	i = good_modulus(i, w);
	j = good_modulus(j, h);
	return getsample_abort_double(x, w, h, pd, i, j, l);
}

// extrapolate by constant (set by calling it with zero sizes)
inline static
double getsample_constant_double(double *x, int w, int h, int pd, int i, int j, int l)
{
	static double value = 0;
	if (w == 0 && h == 0)
		value = *x;
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return value;
	return x[(i+j*w)*pd + l];
}



// test for inclusion of stdlib.h and string.h
#if defined(EXIT_SUCCESS)
#if defined(_STRING_H) || defined(_STRING_H_)
static getsample_operator get_sample_operator(getsample_operator o)
{
	char *option = getenv("GETPIXEL"), *endptr;
	if (!option) return o;
#ifdef NAN
	if (0 == strcmp(option, "nan"      )) return getsample_nan;
#endif//NAN
	if (0 == strcmp(option, "segfault" )) return getsample_error;
	if (0 == strcmp(option, "error"    )) return getsample_error;
	if (0 == strcmp(option, "abort"    )) return getsample_abort;
	if (0 == strcmp(option, "exit"    ))  return getsample_exit;
	if (0 == strcmp(option, "periodic" )) return getsample_per;
	if (0 == strcmp(option, "constant" )) return getsample_0;
	if (0 == strcmp(option, "reflex"   )) return getsample_2;
	if (0 == strcmp(option, "symmetric")) return getsample_2;
	float value = strtof(option, &endptr);
	if (endptr != option) {
		getsample_constant(&value, 0, 0, 0, 0, 0, 0);
		return getsample_constant;
	}
	return getsample_0;
}
#endif//_STRING_H
#endif//EXIT_SUCCESS


inline
static void setsample_0(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		return;
	x[(i+j*w)*pd + l] = v;
}

inline
static void setsample_segf(float *x, int w, int h, int pd, int i, int j, int l,
		float v)
{
	if (i < 0 || i >= w || j < 0 || j >= h || l < 0 || l >= pd)
		*(volatile int*)0 = 0;
	x[(i+j*w)*pd + l] = v;
}

#endif//_GETPIXEL_C2
