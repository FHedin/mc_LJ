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

#define MV_ACC 1
#define MV_REJ -1

int make_MC_moves(ATOM at[], DATA *dat, double *ener)
{
    int i,j,k,accParam=0;
    int acc=0, acc2=0;
    ATOM *at_new;
    int candidate=-1;
    int n_moving=1;
    int *ismoving;
    int mv_direction = -1;
    int progress=dat->nsteps/1000;
    clock_t start,now;

    start=clock();

    at_new=malloc( dat->natom*sizeof *at_new );

    for (i=0; i<(dat->nsteps); i++) //main loop
    {
        for (j=0; j<(dat->natom); j++)
            memcpy(&at_new[j],&at[j],sizeof(ATOM));

        n_moving = 1;
//         n_moving=(int) dat->natom*get_next(dat) + 1;
        ismoving=calloc(n_moving,sizeof *ismoving);

        j=0;
        do
        {
            int redundant=0;
            candidate = (int) dat->natom*get_next(dat);
            for(k=0; k<j; k++)
            {
                if (ismoving[k]==candidate)
                {
                    redundant=1;
                    //                     fprintf(stdout,"Atom %d is redundant.\n",candidate);
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

//                 fprintf(stdout,"%d atoms moving for step %d \n",n_moving,i);
        //         for (j=0;j<n_moving;j++)
        //             fprintf(stdout,"%d\t",ismoving[j]);
        //         fprintf(stdout,"\n");

        double randvec[3] = {0.0}  ;
//         mv_direction = (int)3*get_next(dat);


        k=0;
        do
        {
            j=ismoving[k];
            //             fprintf(stdout,"Moving atom %d\t",j);
            get_vector(dat,mv_direction,randvec);
            at_new[j].x += (dat->d_max)*randvec[0] ;
            at_new[j].y += (dat->d_max)*randvec[1] ;
            at_new[j].z += (dat->d_max)*randvec[2] ;
            k++;
        }
        while(k<n_moving);

        candidate=ismoving[0];
        accParam=apply_Metrop(at,at_new,dat,&candidate,ener,&i);

        if (accParam == MV_ACC)
        {
            acc++;
            acc2++;

            k=0;
            do
            {
                j=ismoving[k];
                memcpy(&at[j],&at_new[j],sizeof(ATOM));
                k++;
            }
            while(k<n_moving);
        }

        adj_dmax(dat,&i,&acc);

//         if ((*ener)/at[0].eps <= dat->E_steepD)
//         {
// 			fprintf(stdout,"Running Steepest Descent at step %d : E = %lf.\n",i,(*ener)/at[0].eps);
//             steepd(at_new,dat);
//             double E_sd = (*get_ENER)(at_new,dat,-1);// /dat->epsilon;
// 			fprintf(stdout,"Steepest Descent done : E = %lf\n",E_sd/at[0].eps);
// 			if ( E_sd/at[0].eps <= dat->E_expected )
//             {
// 				fprintf(stdout,"Best minimum found after %d steps : E = %lf \n",i,E_sd/at[0].eps);
//                 *ener=E_sd;
//                 (*write_traj)(at,dat,i);
//                 return acc2;
//             }
//             else fprintf(stdout,"Minimum not reached : continuing ...\n\n");
//         }

        if (i!=0 && i%io.trsave==0)
            (*write_traj)(at,dat,i);

        //energy check
// 		if (i!=0 && i%1000==0)
// 		{
// 			static double cum_err = 0.0;
// 			double de = (*ener)-((*get_ENER)(at,dat,-1));
// 			cum_err += de;
// 			fprintf(stderr,"At step %d DeltaE is : %g cumulated : %g\n",i,de,cum_err);
// 		}

        if(!is_stdout_redirected && i%progress==0)
        {
            double elapsed=0.0 , remaining =0.0 , pr=0.0 ;
            now = clock();
            elapsed = (double)(now-start)/CLOCKS_PER_SEC;
            pr = 0.1+i/(double)progress/10.0;
            remaining = ((double)progress*elapsed*10.0/(double)i)*(100.0-pr);
            fprintf(stdout,"Progress : %5.1lf %% done \t Estimated remaining time : %8.1lf seconds\r",pr,remaining);
            fflush(stdout);
        }

        free(ismoving);
    }//end of main loop

    free(at_new) ;

    (*write_traj)(at,dat,i);

    return acc2;
}

int apply_Metrop(ATOM at[], ATOM at_new[], DATA *dat, int *candidate, double *ener, int *step)
{
    //return 1 if move accepted, -1 if rejected
    double Eold, Enew, Ediff = 0.;
    double alpha = 0;
    double rejParam = 0;

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
        Eold=(*get_ENER)(at,dat,*candidate);
        Enew=(*get_ENER)(at_new,dat,*candidate);
#ifdef _OPENMP
    }
#endif

    Ediff = Enew - Eold;
    //     fprintf(stdout,"Ediff : %lf\n",Ediff);

    if (Ediff < 0.0)
    {
        *ener+=Ediff;
        if((*step)!=0 && (*step)%io.esave==0)
            fwrite(ener,sizeof(double),1,efile);
        return MV_ACC ;
    }
    else
    {
        rejParam = exp(-dat->beta*Ediff);
        alpha = get_next(dat);
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
