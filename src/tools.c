/**
 * \file tools.c
 *
 * \brief File containing tool functions called across the program
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
#include <string.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "rand.h"
#include "ener.h"
#include "logger.h"

/**
 * \def CONFLICT -1
 * \brief A macro that indicates that there is a steric clash (2 atoms too closed) when the system is randomly generated
 */
#define CONFLICT -1

/**
 * \def NO_CONFLICT 0
 * \brief A macro that indicates that there is no steric clash 
 */
#define NO_CONFLICT 0

/**
 * @brief Fills partially or fully a X,Y,Z vector whith random numbers distributed 
 * in the (-0.5;0.5) range (Reversible Markov Chain, Detailed Balance).
 * 
 * @param dat Common data for simulation
 * @param mv_direction -1 if we want 3 random numbers stored in vec, 0 or 1 or 2 if we want a random X or Y or Z, respectively
 * @param vec A X,Y,Z coordinates vector
 */
void get_vector(DATA *dat, int32_t mv_direction, double vec[3])
{
    vec[0] = 0.0 ;
    vec[1] = 0.0 ;
    vec[2] = 0.0 ;

    if (mv_direction==-1)
    {
        vec[0] = 2.*get_next(dat)-1.;
        vec[1] = 2.*get_next(dat)-1.;
        vec[2] = 2.*get_next(dat)-1.;
    }
    else
        vec[mv_direction] = 2.*get_next(dat)-1.;
    
//     fprintf(stdout,"get_vector : mv_direction is %d\n",mv_direction);
//     fflush(stdout);
    
}

/**
 * @brief Generates an initial cluster of atoms for starting simulation
 * 
 * @param at Atom array
 * @param dat Common data for simulation
 * @param from first atom of the list on which to work
 * @param to last atom of the list on which to work
 * @param mode Building mode : -1 sets coordinates to 9999.9 (infinity), 0 sets all atom at the origin, 1 sets atom at a random position (with constraints)
 */
void build_cluster(ATOM at[], DATA *dat, uint32_t from, uint32_t to, int32_t mode)
{
    uint32_t i = 0 ;
    
//     fprintf(stdout,"Entered in build cluster : from %d to %d mode %d\n",from,to,mode);
//     fflush(stdout);
                    
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
//        double rtot=0.0;
//
//        for (i=0; i<dat->natom; i++)
//            rtot += at[i].ljp.sig;
//        
//        rtot /= (double)dat->natom;
//        rtot = (double)dat->natom*4.0*3.14159*X3(rtot)/3.0;
//        dat->inid = pow(3.0*rtot/(4.0*3.14159),0.333);
        dat->inid = dat->natom/8.0;
        
        for (i=from; i<to; i++)
        {
            
//             fprintf(stdout,"loop %d\n",i);
//             fflush(stdout);
            
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


/**
 * @brief Checks if an atom randomly placed is not far enough from the other ones.
 *          Positioning is valid if the distance is larger than their sigma_i_j LJ interaction term
 * @param at Atom array
 * @param i Atomic index
 * @return NO_CONFLICT if no steric clash, CONFLICT otherwise
 */
int32_t  no_conflict(ATOM at[],uint32_t i)
{
    uint32_t j=0;
    double d=0.0;
    
//     fprintf(stdout,"In no_conflict for %d\n",i);
//     fflush(stdout);
    
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

/**
 * @brief Adjusts the d_max parameter i.e. the maximum move in Angstroems applied to a coordinate
 *      Value is adjusted for reaching d_max_tgt
 *      
 * @param dat Common data for simulation
 * @param step Step at which adjustment is done
 * @param acc The current number of accepted move since last call to this function
 */
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

/**
 * @brief Get the center of mass of the system.
 * As we don't really have mass here this is in fact the barycentre of the system
 * 
 * @param at Atom list
 * @param dat Common data
 * @return The center of mass of the system
 */
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

/**
 * Substracts the center of mass of the system from the coordinates so that 
 * the system is centred at the origin
 * 
 * @param at Atom list
 * @param dat Common data
 */
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
