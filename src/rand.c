/*
 * Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "global.h"
#include "rand.h"
#include "logger.h"

double get_next(DATA *dat)
{
    if (dat->nrn==2048)
    {
#ifdef STDRAND
        uint32_t i;
        for (i=0; i<2048; i++)
            dat->rn[i]=rand()/(double)RAND_MAX;
#else
        dsfmt_fill_array_close_open(&dat->dsfmt,dat->rn,(int32_t)dat->nrn);
#endif
        dat->nrn=0;
    }
    dat->nrn += 1;
    return dat->rn[dat->nrn-1] ;
}

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

// This is a test function for evaluating the quality of the normal random numbers generator.
// generate 'n' normal distributed rand numbers centred around 0 and with spdat->weps as stddev
// numbers are saved in a file norm.dat
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

