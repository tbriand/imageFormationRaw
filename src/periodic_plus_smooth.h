// This program is free software: you can use, modify and/or redistribute it
// under the terms of the simplified BSD License. You should have received a
// copy of this license along this program. If not, see
// <http://www.opensource.org/licenses/bsd-license.html>.
//
// Copyright (C) 2019, Thibaud Briand <briand.thibaud@gmail.com>
// Copyright (C) 2019, Jonathan Vacher <jonathan.vacher@einstein.yu.edu>
// All rights reserved.

/* periodic_component.c contains the functions for computing the periodic
 * plus smooth decomposition of an image */

#ifndef PERIODIC_PLUS_SMOOTH_H
#define PERIODIC_PLUS_SMOOTH_H

void periodic_plus_smooth_decomposition(double *periodic, double *smooth,
                                        const double *in, int w, int h, int pd, double zoom);

#endif
