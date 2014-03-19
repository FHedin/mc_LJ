/*
 * Copyright (c) 2014, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 *
 */

#ifndef MCCLASSIC_H_INCLUDED
#define MCCLASSIC_H_INCLUDED

uint64_t make_MC_moves(ATOM at[], DATA *dat, double *ener);
int32_t apply_Metrop(ATOM at[], ATOM at_new[], DATA *dat, int32_t *candidate, double *ener, uint64_t *step);

#endif // MCCLASSIC_H_INCLUDED
