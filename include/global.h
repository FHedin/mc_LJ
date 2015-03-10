/**
 * \file global.h
 *
 * \brief Header file included everywhere and containing variables / macros / structures of common use
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

#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

#include <stdint.h>
#include <inttypes.h>

#ifndef STDRAND
#include "dSFMT.h"
#endif

#ifndef FILENAME_MAX
#define FILENAME_MAX    4096
#endif

// macro for powers
/**
 * \def X2(a)
 * \brief A macro that squares a number
 */
#define X2(a)   (a)*(a)
/**
 * \def X3(a)
 * \brief A macro that provides the power of 3 of a number
 */
#define X3(a)   X2(a)*(a)
/**
 * \def X4(a)
 * \brief A macro that provides the power of 4 of a number
 */
#define X4(a)   X2(a)*X2(a)
/**
 * \def X6(a)
 * \brief A macro that provides the power of 6 of a number
 */
#define X6(a)   X4(a)*X2(a)
/**
 * \def X12(a)
 * \brief A macro that provides the power of 12 of a number
 */
#define X12(a)  X6(a)*X6(a)

//define where is the null file
#ifdef __unix__
#define NULLFILE "/dev/null"
#else //for MS windows
#define NULLFILE "nul"
#endif

/// if stdout has been redirected to a file from command line call ( -o option)
extern uint32_t is_stdout_redirected;

/// if we decide to use CHARMM units instead of reduced units (not recommended!)
/// 0 = reduced units
/// 1 = CHARMM units
extern uint32_t charmm_units;

#ifdef _OPENMP
extern uint32_t ncpus;
extern uint32_t nthreads;
#endif

/**
 * @brief This structure holds useful variables used across the simulations,
 * it is almost always transmitted from one function to another one .
 */
typedef struct
{
    uint32_t natom ;    ///< Number of atoms
    char method[32];    ///< MC Method string : 'METROP' or 'SPAV' (case insensitive)
    uint64_t nsteps ;   ///< Number of steps as a 64 bits integer to allow really long simulations (i.e. more than 2 billions)

    double d_max ;          ///< Maximum possible distance in Angstroems for a MC move
    uint32_t d_max_when;    ///< when to update d_max (a number of steps)
    double d_max_tgt;       ///< d_max will be tuned in order to reach d_max_tgt percents of moves acceptance

    double inid ;       ///< An initial distance term used when randomly assigning coordinates to atoms when generating a cluster
    double T ;          ///< Temperature : in reduced units, or kcal/mol if charmm_units is 1
    double E_constr;    ///< An energy constraints for avoiding "cluster evaporation" : avoids that atoms go to far from each other : see ener.c for details
    double E_steepD;    ///< A threshold at which starting Steepest Descent minimisation
    double E_expected;  ///< Expected best energy minima of the cluster currently studied ; usually taken from http://www-wales.ch.cam.ac.uk/CCD.html
    double beta;        ///< the inverse temperature used in acceptance criterion

#ifndef STDRAND
    dsfmt_t dsfmt;      ///< A structure used by the dSFMT random numbers generator
    uint32_t *seeds;    ///< An array of seeds used for intialising the dSFMT random numbers generator
#endif
    uint32_t nrn ;      ///< a counter to know how many random numbers from the rn array we have used
    double *rn ;        ///< to avoid calling too often dSFMT, numbers are "cached" i.e. stored in an array ; see rand.c and rand.h
} DATA;

/**
 * @brief This structure holds useful variables used when Spatial Averaging simulations are performed
 * See http://pubs.acs.org/doi/abs/10.1021/ct500529w
 */
typedef struct
{
    uint32_t meps;  ///< Number of sets M_epsilon
    uint32_t neps;  ///< Number of replicated structures N_epsilon
    double weps;    ///< standards deviation of the gaussian distribution
    double *normalNumbs;    ///< to avoid calling too often dSFMT, numbers are "cached" i.e. stored in an array ; see rand.c and rand.h
    uint32_t normalSize;    ///< a counter to know how many random numbers from the normalNumbs array we have used
} SPDAT;

/**
 * @brief A structure holding the Lennard-Jones parameters
 * See http://www.sklogwiki.org/SklogWiki/index.php/Lennard-Jones_model
 */
typedef struct
{
    char sym[4];    ///< atomic symbol
    double sig ;    ///< L-J sigma parameter
    double eps ;    ///< L-J epsilon parameter
} LJPARAMS;

/**
 * @brief A structure representing an atom
 */
typedef struct
{
    double x;   ///< X coordinate
    double y;   ///< Y coordinate
    double z;   ///< Z coordinate
    char sym[4] ;   ///< atomic symbol
    LJPARAMS ljp;    ///< substructure containing LJ parameters
} ATOM;

/**
 * @brief A structure representing the center of mass of a system, simply a point in
 * a 3-Dim space
 */
typedef struct
{
    double cx,cy,cz;
} CM;

#endif // GLOBAL_H_INCLUDED
