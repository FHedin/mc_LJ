/**
 * \file MCclassic.h
 *
 * \brief Header file for MCclassic.c
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

#ifndef MCCLASSIC_H_INCLUDED
#define MCCLASSIC_H_INCLUDED

uint64_t make_MC_moves(ATOM at[], DATA *dat, double *ener);
int32_t apply_Metrop(ATOM at[], ATOM at_new[], DATA *dat, int32_t *candidate, double *ener, uint64_t *step);

#endif // MCCLASSIC_H_INCLUDED
