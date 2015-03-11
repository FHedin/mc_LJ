/**
 * \file minim.h
 *
 * \brief Header file for minim.c
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2011-2015, Florent Hédin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#ifndef MINIM_H_INCLUDED
#define MINIM_H_INCLUDED

void alloc_minim(DATA *dat);
void dealloc_minim();

void steepd(ATOM at[],DATA *dat);
// void steepd_ini(ATOM at[],DATA *dat);
void adjust_alpha(const uint32_t natom, const double grad_old[], const double grad_new[], double *alpha);

#endif // MINIM_H_INCLUDED
