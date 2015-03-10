/**
 * \file rand.h
 *
 * \brief Header file of rand.c file
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

#ifndef RAND_H_INCLUDED
#define RAND_H_INCLUDED

/// get a uniformly distributed random number 
double get_next(DATA *dat);

/// get a normally distributed random number
double get_BoxMuller(DATA *dat, SPDAT *spdat);

/// if we want to test the random numbers generators
void test_norm_distrib(DATA *dat, SPDAT *spdat, uint32_t n);

#endif // RAND_H_INCLUDED
