/*
 * Copyright (c) 2014, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 *
 */

#ifndef MCSPAV_H_INCLUDED
#define MCSPAV_H_INCLUDED

uint64_t launch_SPAV(ATOM at[], DATA *dat, SPDAT *spdat, double *ener);
int32_t apply_SPAV_Criterion(DATA *dat, SPDAT *spdat, ATOM at[], ATOM at_new[],
                             ATOM ***iniArray, ATOM ***finArray, int32_t *candidate,
                             double *ener, uint64_t *currStep);

#endif // MCSPAV_H_INCLUDED
