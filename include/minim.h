/*
 * Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 *
 */

#ifndef MINIM_H_INCLUDED
#define MINIM_H_INCLUDED

void alloc_minim(DATA *dat);
void dealloc_minim();

void steepd(ATOM at[],DATA *dat);
// void steepd_ini(ATOM at[],DATA *dat);
void adjust_alpha(const int natom, const double grad_old[], const double grad_new[], double *alpha);

#endif // MINIM_H_INCLUDED
