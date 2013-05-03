#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "global.h"
#include "MCspav.h"
#include "tools.h"
#include "rand.h"
#include "memory.h"
#include "ener.h"
#include "io.h"

int launch_SPAV(ATOM at[], DATA *dat, SPDAT *spdat, double *ener)
{
    int acc = 0, acc2=0 ;
    int st = 0 ;
    int i,j,k = 0 ;
    int candidate=-1;
    //     int n_moving=1;
    double randvec[3] = {0.};
    int is_accepted = 0 ;
    int mv_direction = -1;
    int progress=dat->nsteps/1000;
    clock_t start,now;

    ATOM *at_new=NULL;
    at_new = calloc(dat->natom,sizeof *at_new);

    ATOM ***iniArray=(ATOM***)calloc_3D(spdat->meps,spdat->neps,dat->natom,sizeof ***iniArray);
    ATOM ***finArray=(ATOM***)calloc_3D(spdat->meps,spdat->neps,dat->natom,sizeof ***finArray);

    for (st=0; st<dat->nsteps; st++)
    {
        for (i=0; i<spdat->meps; i++)
            for (j=0; j<spdat->neps; j++)
                for (k=0; k<dat->natom; k++)
                {
                    memcpy(&at_new[k],&at[k],sizeof(ATOM));
                    memcpy(&iniArray[i][j][k],&at[k],sizeof(ATOM));
                    memcpy(&finArray[i][j][k],&at[k],sizeof(ATOM));
                }

        candidate = (int)dat->natom*get_next(dat);
//         mv_direction = (int)3*get_next(dat);
//         fprintf(stfile,"%d\n",mv_direction);

        k=candidate;

        for (i=0; i<spdat->meps; i++)
            for (j=0; j<spdat->neps; j++)
                //                for (k=0; k<dat->natom; k++)
            {
//                 if (mv_direction==-1 || mv_direction==0) iniArray[i][j][k].x += get_BoxMuller(dat,spdat);
//                 if (mv_direction==-1 || mv_direction==1) iniArray[i][j][k].y += get_BoxMuller(dat,spdat);
//                 if (mv_direction==-1 || mv_direction==2) iniArray[i][j][k].z += get_BoxMuller(dat,spdat);
                iniArray[i][j][k].x += get_BoxMuller(dat,spdat);
                iniArray[i][j][k].y += get_BoxMuller(dat,spdat);
                iniArray[i][j][k].z += get_BoxMuller(dat,spdat);
                memcpy(&finArray[i][j][k],&iniArray[i][j][k],sizeof(ATOM));
            }

        //        for (k=0; k<dat->natom; k++)
        //        {
        randvec[0] = 0.0 ;
        randvec[1] = 0.0 ;
        randvec[2] = 0.0 ;
        get_vector(dat,mv_direction,randvec);
        at_new[k].x += (dat->d_max)*randvec[0] ;
        at_new[k].y += (dat->d_max)*randvec[1] ;
        at_new[k].z += (dat->d_max)*randvec[2] ;
        for (i=0; i<spdat->meps; i++)
            for (j=0; j<spdat->neps; j++)
            {
                finArray[i][j][k].x += (dat->d_max)*randvec[0] ;
                finArray[i][j][k].y += (dat->d_max)*randvec[1] ;
                finArray[i][j][k].z += (dat->d_max)*randvec[2] ;
            }
        //        }

        is_accepted = apply_SPAV_Criterion(dat,spdat,&sts,at,at_new,iniArray,finArray,&candidate,ener,&st);

        if (is_accepted == MV_ACC)
        {
            j=candidate;
            //             fprintf(stdout,"accepted\n");
            acc++;
            acc2++;
            //            for (j=0; j<(dat->natom); j++)
            //            {
            at[j].x = at_new[j].x ;
            at[j].y = at_new[j].y ;
            at[j].z = at_new[j].z ;
            //            }
        }

        adj_dmax(dat,&st,&acc);

//         if (*ener <= dat->E_steepD)
//         {
//             fprintf(stdout,"Running Steepest Descent at step %d .\n",st);
//             steepd(at_new,dat);
//             double E_sd = (*get_ENER)(at_new,dat,-1); // /dat->epsilon;
//             fprintf(stdout,"Steepest Descent done : E = %lf\n",E_sd);
//             if ( E_sd <= dat->E_expected )
//             {
//                 steepd(at_new,dat);
//                 double E_sd = (*get_ENER)(at_new,dat,-1); // /dat->epsilon;
//                 fprintf(stdout,"Best minimum found after %d steps : E = %lf \n",st,E_sd);
//                 *ener=E_sd;
//                 (*write_traj)(at,dat,traj,i);
//                 fprintf(stdout,"STATS : AVERAGES %d :\tdelta : %lf\tsigma : %lf\trejp : %lf\n" , sts.iter , sts.delta_sum/sts.iter , sts.sigma_sum/sts.iter , (sts.delta_sum+sts.sigma_sum/2.0)/sts.iter);
//                 return acc2;
//             }
//             else fprintf(stdout,"Minimum not reached : continuing ...\n\n");
//         }

        if (st!=0 && st%io.trsave==0)
            (*write_traj)(at,dat,st);

        if(!is_stdout_redirected && st%progress==0)
        {
            double elapsed=0.0 , remaining =0.0 , pr=0.0 ;
            now = clock();
            elapsed = (double)(now-start)/CLOCKS_PER_SEC;
            pr = 0.1+st/(double)progress/10.0;
            remaining = ((double)progress*elapsed*10.0/(double)st)*(100.0-pr);
            fprintf(stdout,"Progress : %5.1lf %% done \t Estimated remaining time : %8.1lf seconds\r",pr,remaining);
            fflush(stdout);
        }

    } //END OF MAIN FOR

    // END OF FUNCTION
    free_3D(spdat->meps,spdat->neps,iniArray,finArray,NULL);

    free(at_new);

    (*write_traj)(at,dat,st);

    return acc2;
}

int apply_SPAV_Criterion(DATA *dat, SPDAT *spdat, ATOM at[], ATOM at_new[],
                         ATOM ***iniArray, ATOM ***finArray, int *candidate, double *ener, int *currStep)
{
    int i=0;
    double Eold=0.,Enew=0.,Ediff=0.;

    Eold=(*get_ENER)(at,dat,*candidate);
    Enew=(*get_ENER)(at_new,dat,*candidate);

    Ediff = Enew - Eold;

    if (Ediff < 0.0)
    {
        *ener+=Ediff;
        if((*currStep)!=0 && (*currStep)%io.esave==0)
            fwrite(ener,sizeof(double),1,efile);
        return MV_ACC ;
    }
    else
    {
        int j,k = 0 ;
        double delta = 0. ;
        double sigma = 0. ;
        double rejParam = 0. ;
        double alpha = 0. ;

        double *Smnew=NULL;
        double *Smold=NULL;
        double *deltaM=NULL;

        double **EI=NULL;
        double **EF=NULL;

        Smnew=calloc(spdat->meps,sizeof *Smnew);
        Smold=calloc(spdat->meps,sizeof *Smold);
        deltaM=calloc(spdat->meps,sizeof *deltaM);

        EI=(double**)calloc_2D(spdat->meps,spdat->neps,sizeof **EI);
        EF=(double**)calloc_2D(spdat->meps,spdat->neps,sizeof **EF);

#ifdef _OPENMP
#pragma omp parallel default(shared) firstprivate(i,j)
        {
#pragma omp for schedule(dynamic, 2)
#endif
            for (i=0; i<spdat->meps; i++)
            {
                for (j=0; j<spdat->neps; j++)
                {
                    EI[i][j] = (*get_ENER)(iniArray[i][j],dat,*candidate);
                    EF[i][j] = (*get_ENER)(finArray[i][j],dat,*candidate);
                }
            }
#ifdef _OPENMP
        }
#endif
        for (i=0; i<spdat->meps; i++)
        {
            for (j=0; j<spdat->neps; j++)
            {
                EI[i][j] = exp(-dat->beta*EI[i][j]);
                EF[i][j] = exp(-dat->beta*EF[i][j]);
            }
        }

        for (i=0; i<spdat->meps; i++)
        {
            Smnew[i] = 0.;
            Smold[i] = 0.;
            deltaM[i] = 0.;
            for (j=0; j<spdat->neps; j++)
            {
                Smnew[i] += EF[i][j];
                Smold[i] += EI[i][j];
            }
            deltaM[i] = -1.0*log(Smnew[i]/Smold[i]);
            delta += deltaM[i];
        }

        delta *= 1./spdat->meps;
        for (i = 0 ; i < spdat->meps ; i++)
            sigma += X2(deltaM[i]-delta) ;

        sigma *= 1/(spdat->meps*(spdat->meps-1.0)) ;
        rejParam = delta + sigma/2.0 ;
        rejParam = exp(-dat->beta*rejParam);

        alpha = get_next(dat);

        free(Smnew);
        free(Smold);
        free(deltaM);
        free_2D(spdat->meps,EI,EF,NULL);

        if (alpha < rejParam)
        {
            *ener+=Ediff;
            if((*currStep)!=0 && (*currStep)%io.esave==0)
                fwrite(ener,sizeof(double),1,efile);
            return MV_ACC ;
        }
        else
        {
            if((*currStep)!=0 && (*currStep)%io.esave==0)
                fwrite(ener,sizeof(double),1,efile);
            return MV_REJ ;
        }
    }
}
