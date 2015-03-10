/**
 * \file rand.c
 *
 * \brief File containing functions related to random number generators
 *        If -DSTDRAND the standard C generator is used (NOT RECOMMANDED)
 *        If not by default the dSFMT generator is used (high quality, HIGHLY RECOMMENDED)
 *        Code for dSFMT is included in a subdirectory (unmodified, please avoid modifying excepted if you really know what you do as you may induce a severe bias)
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.h"
#include "rand.h"
#include "logger.h"

/**
 * @brief Call this function for obtaining a uniformly distributed random number in the range (0, 1)
 * 
 * @param dat Common simulation data
 * @return A random number uniformly distributed in the range (0, 1)
 */
double get_next(DATA *dat)
{
    // if array empty of fully used re-fill it
    if (dat->nrn==2048)
    {
#ifdef STDRAND
        uint32_t i;
        for (i=0; i<2048; i++)
            dat->rn[i]=rand()/(double)RAND_MAX;
#else
        dsfmt_fill_array_open_open(&dat->dsfmt,dat->rn,(int32_t)dat->nrn);
#endif
        dat->nrn=0;
    }
    // return a number from the array and increase the counter
    dat->nrn += 1;
    return dat->rn[dat->nrn-1] ;
}

/**
 * @brief Returns a double precision number normally distributed around 0,
 *  following a standard deviation taken from spdat->weps
 *  It uses the Box Muller algorithm, see http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
 * 
 * @param dat Common simulation data
 * @param spdat Spatial-Averaging simulation data
 * 
 * @return A random number normally distributed around 0 and following a standard deviation taken from spdat->weps
 */
double get_BoxMuller(DATA *dat, SPDAT *spdat)
{
    //returns a double normal-distributed around 0, according to a std_dev = spdat->weps
    if (spdat->normalSize==2048)
    {
        uint32_t i;
        double u,v,s;
        for (i=0; i<spdat->normalSize; i++)
        {
            do
            {
                u = 2.*get_next(dat)-1.;
                v = 2.*get_next(dat)-1.;
                s = u*u + v*v;
            }
            while (s >= 1);
            spdat->normalNumbs[i]   = u*spdat->weps*sqrt(-2.*log(s)/s);
            spdat->normalNumbs[i+1] = v*spdat->weps*sqrt(-2.*log(s)/s);
            i++;
        }
        spdat->normalSize=0;
    }
    spdat->normalSize += 1 ;
    //fprintf(stderr,"GAUSS: %lf\n",spdat->normalNumbs[spdat->normalSize-1]);
    return spdat->normalNumbs[spdat->normalSize-1];
}

/**
 * @brief This is a test function for evaluating the quality of the normal random numbers generator.
 *          generates n normal distributed rand numbers centred around 0 and with spdat->weps as stddev
 *          numbers are saved in a file norm.dat
 * 
 * @param dat Common simulation data
 * @param spdat Spatial-Averaging simulation data
 * @param n How many random numbers are required
 */
void test_norm_distrib(DATA *dat, SPDAT *spdat, uint32_t n)
{
    uint32_t i;
    double *norm;
    double mean=0.;
    double sd=0.;

    norm=malloc(n*sizeof *norm);

    FILE *normal=NULL;
    normal=fopen("norm.dat","w");

    for(i=0; i<n; i++)
    {
        norm[i] = get_BoxMuller(dat,spdat);
        mean+=norm[i];
        fprintf(normal,"%lf\n",norm[i]);
    }

    mean/=(double)n;

    for(i=0; i<n; i++)
        sd+=(norm[i]-mean)*(norm[i]-mean);

    sd=sqrt(sd/(double)n);

    fprintf(stdout,"testing normal values\n");
    fprintf(stdout,"mean = %lf\tstddev = %lf\n",mean,sd);

    free(norm);
    norm=NULL;

    fclose(normal);

}

