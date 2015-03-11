/**
 * \file minim.c
 *
 * \brief File containing functions used for minimising the energy of the system
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
#include <math.h>
#include <string.h>

#include "global.h"
#include "ener.h"
#include "logger.h"
#include "minim.h"

static double *fx = NULL ;
static double *fy = NULL ;
static double *fz = NULL ;

static double *fxo = NULL ;
static double *fyo = NULL ;
static double *fzo = NULL ;

void alloc_minim(DATA *dat)
{
    fx  = malloc(dat->natom*sizeof *fx);
    fy  = malloc(dat->natom*sizeof *fy);
    fz  = malloc(dat->natom*sizeof *fz);

    fxo  = malloc(dat->natom*sizeof *fxo);
    fyo  = malloc(dat->natom*sizeof *fyo);
    fzo  = malloc(dat->natom*sizeof *fzo);
}

void dealloc_minim()
{
    free(fx);
    free(fy);
    free(fz);
    free(fxo);
    free(fyo);
    free(fzo);
}

void steepd(ATOM at[],DATA *dat)
{
    uint32_t i=0,/*j=0,*/counter=0;

    double alpha[3] = {1.0e-03,1.0e-03,1.0e-03} ;
    const double prec = 1.0e-05 ;
    double e1=0. , e2=0.;
    double diff;

//     ATOM *at2 = NULL ;

//     at2 = malloc(dat->natom*sizeof *at2);

//     for (i=0; i<(dat->natom); i++)
//         memcpy(&at2[i],&at[i],sizeof(ATOM));

    (*get_DV)(at,dat,fx,fy,fz);
    memcpy(fxo,fx,dat->natom*sizeof(double));
    memcpy(fyo,fy,dat->natom*sizeof(double));
    memcpy(fzo,fz,dat->natom*sizeof(double));

    do
    {
        for (i=0; i<(dat->natom); i++)
        {
            at[i].x -= alpha[0]*fx[i];
            at[i].y -= alpha[1]*fy[i];
            at[i].z -= alpha[2]*fz[i];
        }

        (*get_DV)(at,dat,fx,fy,fz);

//         LOG_PRINT(LOG_DEBUG,"SteepD alpha vector old = %lf %lf %lf\n",alpha[0],alpha[1],alpha[2]);
        adjust_alpha(dat->natom,fxo,fx,alpha);
        adjust_alpha(dat->natom,fyo,fy,alpha+1);
        adjust_alpha(dat->natom,fzo,fz,alpha+2);
//         LOG_PRINT(LOG_DEBUG,"SteepD alpha vector new = %lf %lf %lf\n",alpha[0],alpha[1],alpha[2]);

        memcpy(fxo,fx,dat->natom*sizeof(double));
        memcpy(fyo,fy,dat->natom*sizeof(double));
        memcpy(fzo,fz,dat->natom*sizeof(double));

        e1 = (*get_ENER)(at,dat,-1)/at[0].ljp.eps;
        diff = fabs(e1-e2);
        e2=e1;

        counter++;
    }
    while (diff>prec && counter < 5000);

//     for (i=0; i<(dat->natom); i++)
//         memcpy(&at[i],&at2[i],sizeof(ATOM));
//     free(at2);

    LOG_PRINT(LOG_INFO,"Gradient guided steepest descent took %d steps : diff is %lf \n",counter,diff);
}

// void steepd_ini(ATOM at[],DATA *dat)
// {
//     uint32_t i=0,/*j=0,*/counter=0;
//
//     double step = 0.001 ;
//     double prec = 1.0e-04 ;
//     double diff;
//     double *DV = NULL ;
//
//     DV  = malloc(3*dat->natom*sizeof *DV);
//
//     do
//     {
//         (*get_DV)(at,dat,DV);
//         for (i=0; i<(dat->natom); i++)
//         {
//             at[i].x -= step*DV[i];
//             at[i].y -= step*DV[i+dat->natom];
//             at[i].z -= step*DV[i+2*dat->natom];
//         }
//         //        diff = 0.0;
//         //        for (j=0; j<3*(dat->natom); j++)
//         //            diff += DV[j]*DV[j];
//         //        diff=sqrt(diff);
//         //         fprintf(stderr,"diff = %lf\n",diff);
//         diff = fabs( (*get_ENER)(at,dat,-1)/at[0].ljp.eps - dat->E_expected );
//         //         fprintf(stderr,"diff = %lf\n",diff);
//         counter++;
//     }
//     while (diff>prec && counter <= 1000);
//
//     free(DV);
// }

void adjust_alpha(const uint32_t natom, const double grad_old[], const double grad_new[], double *alpha)
{
    double dot_prod=0.0;
    double norm_old=0.0;
    double norm_new=0.0;
    double angle;

    uint32_t i;

    for(i=0; i<natom; i++)
    {
        dot_prod += grad_old[i]*grad_new[i];
        norm_old += X2(grad_old[i]);
        norm_new += X2(grad_new[i]);
    }
    norm_old = sqrt(norm_old);
    norm_new = sqrt(norm_new);

    angle=acos(dot_prod/(norm_old*norm_new));
    angle*=180.0/3.14159265358979323846;

    if (angle > 60.0)
        *alpha *= 0.5;
    else
        *alpha *= 1.05;
}
