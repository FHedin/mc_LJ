/**
 * \file MCspav.h
 *
 * \brief Header file for MCspav.c
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

#ifndef MCSPAV_H_INCLUDED
#define MCSPAV_H_INCLUDED

uint64_t launch_SPAV(ATOM at[], DATA *dat, SPDAT *spdat, double *ener);
int32_t apply_SPAV_Criterion(DATA *dat, SPDAT *spdat, ATOM at[], ATOM at_new[],
                             ATOM ***iniArray, ATOM ***finArray, int32_t *candidate,
                             double *ener, uint64_t *currStep);

void alloc_SAMC(SPDAT *spdat);
void dealloc_SAMC(SPDAT *spdat);

#endif // MCSPAV_H_INCLUDED
