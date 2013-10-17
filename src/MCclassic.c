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
#include <time.h>

#include "global.h"
#include "MCclassic.h"
#include "tools.h"
#include "rand.h"
#include "ener.h"
#include "io.h"
#include "logger.h"

#define MV_ACC 1
#define MV_REJ -1

uint64_t make_MC_moves(ATOM at[], DATA *dat, double *ener)
{
    uint32_t /*i,*/j,k;
    uint64_t st, acc=0, acc2=0;
    int32_t accParam=0;
    
    ATOM *at_new;
    
//     double *dmax;
    int32_t candidate=-1;
    uint32_t n_moving=   1    ;//dat->natom;
    int32_t *ismoving = NULL;
    //int32_t unicMove=    1;//-1;
    int32_t mv_direction= -1;
    
    double randvec[3] = {0.0,0.0,0.0};
    
//     uint64_t progress=dat->nsteps/1000;
//     clock_t start,now;

//     start=clock();

    at_new=malloc( dat->natom*sizeof *at_new );
    
//     dmax=malloc(dat->natom*sizeof *dmax);
//     for (i=0; i<(dat->natom); i++)
//       dmax[i] = dat->d_max;

    for (st=0; st<(dat->nsteps); st++) //main loop
    {
        for (j=0; j<(dat->natom); j++)
            memcpy(&at_new[j],&at[j],sizeof(ATOM));

//         n_moving=(int) dat->natom*get_next(dat) + 1;
        ismoving=calloc(n_moving,sizeof *ismoving);

        j=0;
        do
        {
            uint32_t redundant=0;
            candidate = (int32_t) (dat->natom*get_next(dat));
            for(k=0; k<j; k++)
            {
                if (ismoving[k]==candidate)
                {
                    redundant=1;
                    break;
                }
            }
            if(redundant)
                continue;
            else
            {
                ismoving[j] = candidate;
                j++;
            }
        }
        while(j<n_moving);

//         fprintf(stderr,"%d atoms moving for step %d \n",n_moving,st);
//         for (j=0;j<n_moving;j++)
//             fprintf(stderr,"%d\t",ismoving[j]);
//         fprintf(stderr,"\n");

        k=0;
        do
        {
            j = (uint32_t) ismoving[k];
//            mv_direction = (int)3*get_next(dat);
            get_vector(dat,mv_direction,randvec);
            
            at_new[j].x += (dat->d_max)*randvec[0] ;
            at_new[j].y += (dat->d_max)*randvec[1] ;
            at_new[j].z += (dat->d_max)*randvec[2] ;
            k++;
        }
        while(k<n_moving);
        
        accParam=apply_Metrop(at,at_new,dat,&ismoving[0],ener,&st);
//         accParam=apply_Metrop(at,at_new,dat,&unicMove,ener,&st);

        if (accParam == MV_ACC)
        {
            acc++;
            acc2++;

            k=0;
            do
            {
                j = (uint32_t) ismoving[k];
                memcpy(&at[j],&at_new[j],sizeof(ATOM));
                k++;
            }
            while(k<n_moving);
        }

        if (dat->d_max_when != 0)
          adj_dmax(dat,&st,&acc);

//         if ((*ener)/at[0].ljp.eps <= dat->E_steepD)
//         {
//           fprintf(stdout,"Running Steepest Descent at step %d : E = %lf.\n",st,(*ener)/at[0].ljp.eps);
//           steepd(at_new,dat);
//           double E_sd = (*get_ENER)(at_new,dat,-1);
//           fprintf(stdout,"Steepest Descent done (step %d): E = %lf\n",st,E_sd/at[0].ljp.eps);
//           if ( fabs(E_sd/at[0].ljp.eps - dat->E_expected ) <= 1.0e-04 )
//           {
//             fprintf(stdout,"Best minimum found after %d steps : E = %lf \n",st,E_sd/at[0].ljp.eps);
//             *ener=E_sd;
//             fwrite(ener,sizeof(double),1,efile);
//             (*write_traj)(at_new,dat,st);
//             return acc2;
//           }
//           else 
//               fprintf(stdout,"Minimum not reached : continuing ...\n\n");
//         }

        if (st!=0 && st%io.trsave==0)
        {
            steepd(at_new,dat);
            double E_sd = (*get_ENER)(at_new,dat,-1);
            fprintf(stdout,"Steepest Descent done (step %"PRIu64"): E = %.3lf\n",st,E_sd/at[0].ljp.eps);
            (*write_traj)(at,dat,st);
        }

        //energy check
//         if (st!=0 && st%1000==0)
//         {
//           static double cum_err = 0.0;
//           double de = (*ener)-((*get_ENER)(at,dat,-1));
//           cum_err += de;
// //           fprintf(stderr,"At step %d DeltaE is : %g cumulated : %g\n",st,de,cum_err);
//         }

//        if(!is_stdout_redirected && st%progress==0)
//        {
//            double elapsed=0.0 , remaining =0.0 , pr=0.0 ;
//            int n = 1;
//#ifdef _OPENMP
//            n = nthreads;
//#endif
//
//            now = clock();
//            elapsed = (double)(now-start)/(n*CLOCKS_PER_SEC);
//            pr = 0.1+st/(double)progress/10.0;
//            remaining = ((double)progress*elapsed*10.0/(double)st)*(100.0-pr);
//            fprintf(stdout,"Progress : %5.1lf %% done \t Estimated remaining time : %8.1lf seconds\r",pr,remaining);
//            fflush(stdout);
//        }

        free(ismoving);
        
    }//end of main loop

    free(at_new) ;

    (*write_traj)(at,dat,st);

    return acc2;
}

int32_t apply_Metrop(ATOM at[], ATOM at_new[], DATA *dat, int32_t *candidate, double *ener, uint64_t *step)
{
    //return 1 if move accepted, -1 if rejected
    double Eold=0.0, Enew=0.0, Ediff=0.0;
    double EconstrOld=0.0,EconstrNew=0.0,EconstrDiff=0.0;
    double alpha = 0.;
    double rejParam = 0.;

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
        Eold=(*get_ENER)(at,dat,*candidate);
        EconstrOld=dat->E_constr;

        Enew=(*get_ENER)(at_new,dat,*candidate);
        EconstrNew=dat->E_constr;
#ifdef _OPENMP
    }
#endif

    Ediff = (Enew - Eold) ;
    EconstrDiff = (EconstrNew - EconstrOld) ;
    
//     fprintf(stderr,"Ediff : %lf \t Econstrdiff : %lf \t",Ediff,EconstrDiff);

    if ( (Ediff + EconstrDiff) < 0.0 )
    {
        *ener+=Ediff;
        if((*step)!=0 && (*step)%io.esave==0)
            fwrite(ener,sizeof(double),1,efile);
//         fprintf(stderr,"\n");
        return MV_ACC ;
    }
    else
    {
        rejParam = exp(-dat->beta*(Ediff + EconstrDiff));
        alpha = get_next(dat);
        
//         fprintf(stderr,"alpha : %lf ; \t rejp : %lf\n",alpha,rejParam);
        
        if (alpha < rejParam)
        {
            *ener+=Ediff;
            if((*step)!=0 && (*step)%io.esave==0)
                fwrite(ener,sizeof(double),1,efile);
            return MV_ACC ;
        }
        else
        {
            if((*step)!=0 && (*step)%io.esave==0)
                fwrite(ener,sizeof(double),1,efile);
            return MV_REJ ;
        }
    }
}
