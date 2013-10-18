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
#include <string.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "rand.h"
#include "ener.h"
#include "logger.h"

#define CONFLICT -1
#define NO_CONFLICT 0

void get_vector(DATA *dat, int32_t mv_direction, double vec[3])
{
    vec[0] = 0.0 ;  vec[1] = 0.0 ;  vec[2] = 0.0 ;
    
    if (mv_direction==-1)
    {
        vec[0] = 2.*get_next(dat)-1.;
        vec[1] = 2.*get_next(dat)-1.;
        vec[2] = 2.*get_next(dat)-1.;
    }
    else
        vec[mv_direction] = 2.*get_next(dat)-1.;
}

/*
 * Generates a cluster for the atoms
 */
void build_cluster(ATOM at[], DATA *dat, uint32_t from, uint32_t to, int32_t mode)
{
    uint32_t i = 0 ;
    if (mode==-1)	//infinite initialisation mode
    {
        for (i=from; i<to; i++)
            at[i].x=at[i].y=at[i].z=9999.9;
    }
    else if (mode==0)	//zero everywhere : all at origin
    {
        for (i=from; i<to; i++)
            at[i].x=at[i].y=at[i].z=0.0;
    }
    else if (mode==1)	//random mode
    {
        double randvec[3] = {0.0} ;
        double rtot=0.0;

        for (i=0; i<dat->natom; i++)
            rtot += at[i].ljp.sig;
        rtot /= (double)dat->natom;
        rtot = (double)dat->natom*4.0*3.14159*X3(rtot)/3.0;
        dat->inid = pow(3.0*rtot/(4.0*3.14159),0.333);

        for (i=from; i<to; i++)
        {
            do
            {
                get_vector(dat,-1,randvec);
                at[i].x = dat->inid*randvec[0];
                at[i].y = dat->inid*randvec[1];
                at[i].z = dat->inid*randvec[2];
            }
            while( (X2(at[i].x)+X2(at[i].y)+X2(at[i].z))>X2(dat->inid) || no_conflict(at,i) != NO_CONFLICT );
        }
    }
}

//returns 0 if no steric clash, -1 otherwise
int32_t  no_conflict(ATOM at[],uint32_t i)
{

    uint32_t j=0;
    double d=0.0;
    for (j=0; j<i; j++)
    {
        d = X2(at[i].x-at[j].x) +  X2(at[i].y-at[j].y) + X2(at[i].z-at[j].z) ;
        d = sqrt(d);
        if (d < ((at[i].ljp.sig+at[j].ljp.sig)/2.0))
        {
//             fprintf(stderr,"[Info] Atoms %3d and %3d too close for starting configuration : generating new coordinates for atom %3d\n",j,i,i);
            LOG_PRINT(LOG_INFO,"Atoms %d and %d too close for starting configuration : generating new coordinates for atom %3d\n",j,i,i);
            return CONFLICT;
        }
    }
    return NO_CONFLICT;
}

void adj_dmax(DATA *dat, uint64_t *step, uint64_t *acc)
{
    if (*step != 0 && *step%dat->d_max_when==0)
    {
        double ratio = (double)*acc/(double)dat->d_max_when;
        double new_dmax = dat->d_max;

        if (ratio > dat->d_max_tgt/100)
            new_dmax *= 1.10 ;
        else
            new_dmax *= 0.90 ;
        *acc = 0 ;
        if (new_dmax > 1.0) new_dmax = 1.0;
        if (new_dmax < 0.01) new_dmax = 0.01;
        
        LOG_PRINT(LOG_INFO,"d_max update at step %"PRIu64" : ratio = %lf ; old_dmax = %lf ; new_dmax = %lf\n",*step,ratio,dat->d_max,new_dmax);
        dat->d_max = new_dmax;
    }
}

CM getCM(ATOM at[],DATA *dat)
{
	CM cm;
	cm.cx=0.0;
	cm.cy=0.0;
	cm.cz=0.0;

	for(uint32_t i=0; i<dat->natom; i++)
	{
		cm.cx += at[i].x;
		cm.cy += at[i].y;
		cm.cz += at[i].z;
	}

	cm.cx /= dat->natom;
    	cm.cy /= dat->natom;
    	cm.cz /= dat->natom;

	return cm;
}

void recentre(ATOM at[],DATA *dat)
{
    CM cm = getCM(at,dat);

    for(uint32_t i=0; i<dat->natom; i++)
    {
        at[i].x -= cm.cx;
        at[i].y -= cm.cy;
        at[i].z -= cm.cz;
    }
}
