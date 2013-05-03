#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "rand.h"
#include "ener.h"

#define CONFLICT -1
#define NO_CONFLICT 0

void get_vector(DATA *dat, int mv_direction, double vec[3])
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
void build_cluster(ATOM at[], DATA *dat, int from, int to, int mode)
{
    int i = 0 ;
    if (mode==-1)	//initialisation mode
    {
        for (i=from; i<to; i++)
            at[i].x=at[i].y=at[i].z=9999.9;
    }
    else if (mode==0)	//zero everywhere
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
int  no_conflict(ATOM at[],int i)
{

    int j=0;
    double d=0.0;
    for (j=0; j<i; j++)
    {
        d = X2(at[i].x-at[j].x) +  X2(at[i].y-at[j].y) + X2(at[i].z-at[j].z) ;
        d = sqrt(d);
        if (d < ((at[i].ljp.sig+at[j].ljp.sig)/2.0))
        {
            fprintf(stderr,"Atoms %3d and %3d too close for starting configuration : generating new coordinates for atom %3d\n",j,i,i);
            return CONFLICT;
        }
    }
    return NO_CONFLICT;
}

void steepd(ATOM at[],DATA *dat)
{
    int i=0,j=0,counter=0;

    double step = 0.001 ;
    double prec = 1.0e-04 ;
    double diff;
    double *DV = NULL ;
    ATOM *at2 = NULL ;

    DV  = malloc(3*dat->natom*sizeof *DV);
    at2 = malloc(dat->natom*sizeof *at2);

    for (i=0; i<(dat->natom); i++)
        memcpy(&at2[i],&at[i],sizeof(ATOM));

    do
    {
        (*get_DV)(at2,dat,DV);
        for (i=0; i<(dat->natom); i++)
        {
            at2[i].x -= step*DV[i];
            at2[i].y -= step*DV[i+dat->natom];
            at2[i].z -= step*DV[i+2*dat->natom];
        }
//        diff = 0.0;
//        for (j=0; j<3*(dat->natom); j++)
//            diff += DV[j]*DV[j];
//        diff=sqrt(diff);
//         fprintf(stderr,"diff = %lf\n",diff);
        diff = fabs( (*get_ENER)(at2,dat,-1)/at2[0].ljp.eps - dat->E_expected );
        fprintf(stderr,"diff = %lf\n",diff);
        counter++;
    }
    while (diff>prec && counter <= 10000);
    
    for (i=0; i<(dat->natom); i++)
        memcpy(&at[i],&at2[i],sizeof(ATOM));

    free(DV);
    free(at2);
}

void adj_dmax(DATA *dat, int *step, int *acc)
{
    if (*step != 0 && *step%dat->d_max_when==0)
    {
        double ratio = (double)*acc/(double)dat->d_max_when;
//        double ratio = (double)*acc/(double)*step;
        fprintf(stderr,"d_max update at step %d : ratio = %lf ; old = %lf ; ",*step,ratio,dat->d_max);
        if (ratio > dat->d_max_tgt/100)
            dat->d_max *= 1.10 ;
        else
            dat->d_max *= 0.90 ;
        *acc = 0 ;
        if (dat->d_max > 1.0) dat->d_max = 1.0;
        if (dat->d_max < 0.01) dat->d_max = 0.01;
        fprintf(stderr,"new = %lf\n",dat->d_max);
    }
}

void recentre(ATOM at[],DATA *dat)
{
    double cx=0., cy=0., cz=0. ;
    
    for(int i=0; i<dat->natom; i++)
    {
        cx += at[i].x;
        cy += at[i].y;
        cz += at[i].z;
    }
    
    cx /= dat->natom;
    cy /= dat->natom;
    cz /= dat->natom;
    
    for(int i=0; i<dat->natom; i++)
    {
        at[i].x -= cx;
        at[i].y -= cy;
        at[i].z -= cz;
    }
}