/**
 * \file ener.h
 *
 * \brief Header file for ener.c
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

#ifndef ENER_H_INCLUDED
#define ENER_H_INCLUDED

/*
 * Lennard Jones : divide by this for transforming
 * R_min in sigma ( it is 2^(1/6) )
 * */
#define RTOSIG 1.122462048309373

/*
 * Boltzmann cst in charmm units
 */
#define KBCH 1.98719e-03

/*
 * for converting energy
 */
#define JTOCAL      0.239005736     // Joules to Calories
#define CM1TOKJM    1.1963e-02    // 1 cm-1 in kJ/mol

// pointers to the desired energy and force functions
double (*get_ENER)(ATOM at[], DATA *dat, int32_t candidate);
void   (*get_DV)(ATOM at[], DATA *dat, double fx[], double fy[], double fz[]);

// ener and force for lennard-jones
double get_LJ_V(ATOM at[], DATA *dat, int32_t candidate);
void get_LJ_DV(ATOM at[], DATA *dat, double fx[], double fy[], double fz[]);

// ener for aziz potential
double get_AZIZ_V(ATOM at[], DATA *dat, int32_t candidate);

// those 3 functions returns energy in cm-1 !!
double aziz_ne_ne(double r);
double aziz_ar_ne(double r);
double aziz_ar_ar(double r);

// constraint for avoiding cluster evaporation
double getExtraPot(double d2, double sig, double eps);

#endif // ENER_H_INCLUDED
