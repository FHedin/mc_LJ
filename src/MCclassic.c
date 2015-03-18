/**
 * \file MCclassiv.c
 *
 * \brief Functions for running Metropolis Monte Carlo simulations
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2011-2015, Florent Hedin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
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
#include "minim.h"
#include "io.h"
#include "logger.h"

/**
 * \def MV_ACC 1
 * \brief A macro indicating that the mc move was accepted
 */
#define MV_ACC 1
/**
 * \def MV_REJ -1
 * \brief A macro indicating that the mc move was rejected
 */
#define MV_REJ -1

/**
 * @brief This is the core function for Metropolis MC simulation 
 *        where the main loop is located, 
 * 
 * @param at Atom list
 * @param dat Common data
 * @param ener Variable containing total energy of the system
 * 
 * @return The number of moves accepted
 */
uint64_t make_MC_moves(ATOM at[], DATA *dat, double *ener)
{
    uint32_t /*i,*/j,k;
    uint64_t st, acc=0, acc2=0;
    int32_t accParam=0;

    // a copy of the atom list
    ATOM *at_new;

    //the candidate moving atom
    int32_t candidate =-1;
    //number of simultaneously moving atoms
    uint32_t n_moving = 1;
    //a list of moving atoms
    int32_t *ismoving = NULL;
    //move direction : 0=x 1=y 2=z -1=all
    int32_t mv_direction= -1;

    double randvec[3] = {0.0,0.0,0.0};

    at_new=malloc(dat->natom*sizeof *at_new);
    ismoving=calloc(dat->natom,sizeof *ismoving);

    // main iteration over all steps
    for (st=1; st<=(dat->nsteps); st++) //main loop
    {

        LOG_PRINT(LOG_DEBUG,"----------------------"
                  " STEP %"PRIu64" ----------------------\n",st);

        for (j=0; j<(dat->natom); j++)
            memcpy(&at_new[j],&at[j],sizeof(ATOM));

        // choose how many atoms will move at this step
//         n_moving=(int) dat->natom*get_next(dat) + 1;

        // randomly choose atom(s) moving at this step
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

        LOG_PRINT(LOG_DEBUG,"%d atoms moving for step %"PRIu64" --> ",n_moving,st);
        for (j=0; j<n_moving; j++)
            LOG_PRINT_SHORT(LOG_DEBUG,"%d ",ismoving[j]);
        LOG_PRINT_SHORT(LOG_DEBUG,"\n");

        //apply random moves to moving atoms
        k=0;
        do
        {
            j = (uint32_t) ismoving[k];
            //uncomment if anisotropic moves required (i.e. in only one direction)
//            mv_direction = (int)3*get_next(dat);
            get_vector(dat,mv_direction,randvec);

            at_new[j].x += (dat->d_max)*randvec[0] ;
            at_new[j].y += (dat->d_max)*randvec[1] ;
            at_new[j].z += (dat->d_max)*randvec[2] ;
            k++;
        }
        while(k<n_moving);

        //get acceptance criterion
        accParam=apply_Metrop(at,at_new,dat,&ismoving[0],ener,&st);

        //if accepted
        if (accParam == MV_ACC)
        {
            // increase acceptance counters
            acc++;
            acc2++;

            //copy new coordinates of the moving atom(s)
            k=0;
            do
            {
                j = (uint32_t) ismoving[k];
                memcpy(&at[j],&at_new[j],sizeof(ATOM));
                k++;
            }
            while(k<n_moving);
        }

        //if required adjust dmax
        if (dat->d_max_when != 0)
            adj_dmax(dat,&st,&acc);

        // old code should consider to remove
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

        
        //if necessary save coordinates
        //energy is minimised so we get something similar to D Wales Basin Hopping
	uint32_t sddone = 0;
	double E_sd = 0.;
        if (st!=0 && st%io.trsave==0)
        {
            steepd(at_new,dat);
            sddone=1;
            E_sd = (*get_ENER)(at_new,dat,-1);
            fprintf(stdout,"Steepest Descent done (step %"PRIu64"): E = %.3lf\n",st,E_sd);
            //(*write_traj)(at,dat,st);
	    (*write_traj)(at_new,dat,st);
        }
        
        //if necessary save energy and run steepest descent if not done yet
        if (st!=0 && st%io.esave==0)
	{
	    if(!sddone)
	    {
	      steepd(at_new,dat);
	      sddone=1;
	      E_sd = (*get_ENER)(at_new,dat,-1);
	      fprintf(stdout,"Steepest Descent done (step %"PRIu64"): E = %.3lf\n",st,E_sd);
	    }
	    fwrite(&E_sd,sizeof(double),1,efile);
	}

    }//end of main loop

    free(at_new) ;
    free(ismoving);

    (*write_traj)(at,dat,st);

    return acc2;
}

/**
 * @bried This function is in charge of checking the energy difference between the new and old atomic configurations
 *          and then return if the move is accepted or rejected
 * 
 * @param at Atom list
 * @param at_new Modifed atom lsit
 * @param dat Common data
 * @param candidate List of atoms that were moving
 * @param ener Variable where energy difference will be stored
 * @param step The current simulation step
 * 
 * @return MV_ACC or MV_REJ if the move is either accepted or rejected
 */
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

    LOG_PRINT(LOG_DEBUG,"Ediff : %lf \t Econstrdiff : %lf \n",Ediff,EconstrDiff);

    if ( (Ediff + EconstrDiff) < 0.0 )
    {
        *ener+=Ediff;

        LOG_PRINT(LOG_DEBUG,"MOVE ACCEPTED\n");

        return MV_ACC ;
    }
    else
    {
        rejParam = exp(-dat->beta*(Ediff + EconstrDiff));
        alpha = get_next(dat);

        LOG_PRINT(LOG_DEBUG,"alpha : %lf ; \t rejp : %lf\n",alpha,rejParam);

        if (alpha < rejParam)
        {
            *ener+=Ediff;

            LOG_PRINT(LOG_DEBUG,"MOVE ACCEPTED\n");

            return MV_ACC ;
        }
        else
        {

            LOG_PRINT(LOG_DEBUG,"MOVE REJECTED\n");

            return MV_REJ ;
        }
    }
}
