#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

#define X2(a)   (a)*(a)
#define X3(a)   X2(a)*(a)
#define X4(a)   X2(a)*X2(a)
#define X6(a)   X4(a)*X2(a)
#define X12(a)  X6(a)*X6(a)

#define MV_ACC 1
#define MV_REJ -1

#include <stdint.h>
#include <inttypes.h>

#ifndef STDRAND
#include "dSFMT.h"
#endif

extern uint32_t is_stdout_redirected;
extern uint32_t charmm_units;

#ifdef _OPENMP
extern uint32_t ncpus;
extern uint32_t nthreads;
#endif

typedef struct
{
    uint32_t natom ;
    char method[32];
    uint64_t nsteps ;

    double d_max ;
    uint32_t d_max_when;
    double d_max_tgt;

    double inid ;
    double T ;
    double E_constr;
    double E_steepD;
    double E_expected;
    double beta;
    
#ifndef STDRAND
    dsfmt_t dsfmt;
    uint32_t *seeds;
#endif
    uint32_t nrn ;
    double *rn ;
} DATA;

typedef struct
{
    uint32_t meps;
    uint32_t neps;
    double weps;
    double *normalNumbs;
    uint32_t normalSize;
} SPDAT;

typedef struct
{
    char sym[4];
    double sig ;
    double eps ;
} LJPARAMS;

typedef struct
{
    double x,y,z ;
    char sym[4] ;
    LJPARAMS ljp;
} ATOM;

typedef struct
{
	double cx,cy,cz;
} CM;

#endif // GLOBAL_H_INCLUDED
