/*
 * Copyright (c) 2014, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 *
 */

#ifndef RAND_H_INCLUDED
#define RAND_H_INCLUDED

double get_next(DATA *dat);
double get_BoxMuller(DATA *dat, SPDAT *spdat);
void test_norm_distrib(DATA *dat, SPDAT *spdat, uint32_t n);

#endif // RAND_H_INCLUDED
