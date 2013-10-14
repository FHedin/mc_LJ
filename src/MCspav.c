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

uint64_t launch_SPAV(ATOM at[], DATA *dat, SPDAT *spdat, double *ener)
{
    uint64_t acc=0, acc2=0 ;
    uint64_t st=0 ;
    
    uint32_t i=0,j=0,k=0,l=0;
    
    int32_t candidate = -1;
    uint32_t n_moving = 1;
    int32_t *ismoving = NULL;
    //int32_t unicMove = -1;
    int32_t mv_direction = -1;
    
    double randvec[3] = {0.0,0.0,0.0};

    int32_t is_accepted = 0 ;
//     uint64_t progress=dat->nsteps/1000;
//     clock_t start,now;
    
//     start=clock();

    ATOM *at_new=NULL;
    at_new=malloc( dat->natom*sizeof *at_new );

    ATOM ***iniArray=(ATOM***)calloc_3D(spdat->meps,spdat->neps,dat->natom,sizeof ***iniArray);
    ATOM ***finArray=(ATOM***)calloc_3D(spdat->meps,spdat->neps,dat->natom,sizeof ***finArray);

    for (st=0; st<dat->nsteps; st++)
    {
        for (i=0; i<spdat->meps; i++)
        {
            for (j=0; j<spdat->neps; j++)
            {
                for (k=0; k<dat->natom; k++)
                {
                    memcpy(&(at_new[k]),&(at[k]),sizeof(ATOM));
                    memcpy(&(iniArray[i][j][k]),&(at[k]),sizeof(ATOM));
                    memcpy(&(finArray[i][j][k]),&(at[k]),sizeof(ATOM));
                    
//                    at_new[k].x = at[k].x;     
//                    at_new[k].y = at[k].y;  
//                    at_new[k].y = at[k].y;
//                    strcpy(at_new[k].sym,at[k].sym);
//                    at_new[k].ljp.eps = at[k].ljp.eps;
//                    at_new[k].ljp.sig = at[k].ljp.sig;
//                    strcpy(at_new[k].ljp.sym,at[k].ljp.sym);
                }
            }
        }

        //n_moving=(int) dat->natom*get_next(dat) + 1;
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
        
//        fprintf(stderr,"%d atoms moving for step %d \n",n_moving,st);
//        for (j=0;j<n_moving;j++)
//          fprintf(stderr,"%d\t",ismoving[j]);
//        fprintf(stderr,"\n");

        for (l=0; l<n_moving; l++)
        {
          k = (uint32_t) ismoving[l];
//          k=ismoving[0];
//          mv_direction = (int)3*get_next(dat);
          get_vector(dat,mv_direction,randvec);
          
          at_new[k].x += (dat->d_max)*randvec[0] ;
          at_new[k].y += (dat->d_max)*randvec[1] ;
          at_new[k].z += (dat->d_max)*randvec[2] ;
            
          for (i=0; i<spdat->meps; i++)
          {
            for (j=0; j<spdat->neps; j++)
            {
              iniArray[i][j][k].x += get_BoxMuller(dat,spdat);
              iniArray[i][j][k].y += get_BoxMuller(dat,spdat);
              iniArray[i][j][k].z += get_BoxMuller(dat,spdat);
              
//               memcpy(&finArray[i][j][k],&at_new[k],sizeof(ATOM));
              
//              finArray[i][j][k].x = at_new[k].x +  (dat->d_max)*get_BoxMuller(dat,spdat);
//              finArray[i][j][k].y = at_new[k].x +  (dat->d_max)*get_BoxMuller(dat,spdat);
//              finArray[i][j][k].z = at_new[k].x +  (dat->d_max)*get_BoxMuller(dat,spdat);
              finArray[i][j][k].x = iniArray[i][j][k].x + (dat->d_max)*randvec[0];
              finArray[i][j][k].y = iniArray[i][j][k].y + (dat->d_max)*randvec[1];
              finArray[i][j][k].z = iniArray[i][j][k].z + (dat->d_max)*randvec[2];
            }
          } 
          
        }

        is_accepted = apply_SPAV_Criterion(dat,spdat,at,at_new,iniArray,finArray,&ismoving[0],ener,&st);
//         is_accepted = apply_SPAV_Criterion(dat,spdat,at,at_new,iniArray,finArray,&unicMove,ener,&st);

//        for (k=0; k<dat->natom; k++)
//        {
//            fprintf(stderr,"at[%d]\t%lf\t%lf\t%lf\n",k,at[k].x,at[k].y,at[k].z);
//            fprintf(stderr,"at_new[%d]\t%lf\t%lf\t%lf\n",k,at_new[k].x,at_new[k].y,at_new[k].z);
//            for (i=0; i<spdat->meps; i++)
//            {
//                for (j=0; j<spdat->neps; j++)
//                {
//                        fprintf(stderr,"iniArray[%d][%d][%d]\t%lf\t%lf\t%lf\n",i,j,k,iniArray[i][j][k].x,iniArray[i][j][k].y,iniArray[i][j][k].z);
//                        fprintf(stderr,"finArray[%d][%d][%d]\t%lf\t%lf\t%lf\n",i,j,k,finArray[i][j][k].x,finArray[i][j][k].y,finArray[i][j][k].z);
//                }
//            }
//            fprintf(stderr,"\n");
//        }
        
        if (is_accepted == MV_ACC)
        {
            acc++;
            acc2++;
            for (l=0; l<n_moving; l++)
            {
                j = (uint32_t) ismoving[l];
                
                at[j].x = at_new[j].x ;
                at[j].y = at_new[j].y ;
                at[j].z = at_new[j].z ;
            }
        }
        if (dat->d_max_when != 0)
            adj_dmax(dat,&st,&acc);

//         if ((*ener)/at[0].ljp.eps <= dat->E_steepD)
//         {
//           fprintf(stdout,"Running Steepest Descent at step %d : E = %lf.\n",st,(*ener)/at[0].ljp.eps);
//           steepd(at_new,dat);
//           double E_sd = (*get_ENER)(at_new,dat,-1);
//           fprintf(stdout,"Steepest Descent done : E = %lf\n",E_sd/at[0].ljp.eps);
//           if ( fabs(E_sd/at[0].ljp.eps - dat->E_expected) <= 1.0e-04 )
//           {
//             fprintf(stdout,"Best minimum found after %d steps : E = %lf \n",st,E_sd/at[0].ljp.eps);
//             *ener=E_sd;
//             (*write_traj)(at_new,dat,st);
//             return acc2;
//           }
//           else fprintf(stdout,"Minimum not reached : continuing ...\n\n");
//         }

        if (st!=0 && st%io.trsave==0)
        {
            steepd(at_new,dat);
            double E_sd = (*get_ENER)(at_new,dat,-1);
            fprintf(stdout,"Steepest Descent done (step %"PRIu64"): E = %.3lf\n",st,(E_sd/at[0].ljp.eps) );
            (*write_traj)(at,dat,st);
        }
        
        //energy check
        if (st!=0 && st%1000==0)
        {
          static double cum_err = 0.0;
          double de = (*ener)-((*get_ENER)(at,dat,-1));
          cum_err += de;
          fprintf(stderr,"At step %"PRIu64" DeltaE is : %g cumulated : %g\n",st,de,cum_err);
        }

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
        
    } //END OF MAIN FOR

    free_3D(spdat->meps,spdat->neps,iniArray,finArray,NULL);

    free(at_new);

    (*write_traj)(at,dat,st);

    return acc2;
}

int32_t apply_SPAV_Criterion(DATA *dat, SPDAT *spdat, ATOM at[], ATOM at_new[],
                          ATOM ***iniArray, ATOM ***finArray, int32_t *candidate, double *ener, uint64_t *currStep)
{
    double Eold=0.,Enew=0.,Ediff=0.;
    double EconstrOld=0.0,EconstrNew=0.0,EconstrDiff=0.0;

    Eold=(*get_ENER)(at,dat,*candidate);
    EconstrOld=dat->E_constr;
    
    Enew=(*get_ENER)(at_new,dat,*candidate);
    EconstrNew=dat->E_constr;
    
    

    Ediff = (Enew - Eold) ;
    EconstrDiff = (EconstrNew - EconstrOld) ;
    
//     fprintf(stderr,"Ediff : %lf \t Econstrdiff : %lf \t",Ediff,EconstrDiff);

    if ( (Ediff + EconstrDiff) < 0.0 )
    {
        *ener+=Ediff;
        if((*currStep)!=0 && (*currStep)%io.esave==0)
            fwrite(ener,sizeof(double),1,efile);
//         fprintf(stderr,"\n");
        return MV_ACC ;
    }
    else
    {
        uint32_t i,j/*,k*/ = 0 ;
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
//            fputs("\n",stderr);
            for (i=0; i<spdat->meps; i++)
            {
                for (j=0; j<spdat->neps; j++)
                {
                    EI[i][j] =  (*get_ENER)(iniArray[i][j],dat,*candidate);
                    EI[i][j] += dat->E_constr;
                    
                    EF[i][j] = (*get_ENER)(finArray[i][j],dat,*candidate);
                    EF[i][j] += dat->E_constr;
                    
//                    fprintf(stderr,"EI[%d][%d]=%lf \t EF[%d][%d]=%lf \n",i,j,EI[i][j],i,j,EF[i][j]);
                }
            }
//            fputs("\n",stderr);
#ifdef _OPENMP
        }
#endif
//        fputs("\n",stderr);
        for (i=0; i<spdat->meps; i++)
        {
            for (j=0; j<spdat->neps; j++)
            {
                EI[i][j] = exp(-dat->beta*EI[i][j]);
                EF[i][j] = exp(-dat->beta*EF[i][j]);
//                fprintf(stderr,"EI[%d][%d]=%lf \t EF[%d][%d]=%lf \n",i,j,EI[i][j],i,j,EF[i][j]);
            }
        }
//        fputs("\n",stderr);

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
        
//         fprintf(stderr,"alpha : %lf ; \t rejp : %lf",alpha,rejParam);

        free(Smnew);
        free(Smold);
        free(deltaM);
        free_2D(spdat->meps,EI,EF,NULL);

        if (alpha < rejParam)
        {
            *ener+=Ediff;
            if((*currStep)!=0 && (*currStep)%io.esave==0)
                fwrite(ener,sizeof(double),1,efile);
//             fprintf(stderr,"\n");
            return MV_ACC ;
        }
        else
        {
            if((*currStep)!=0 && (*currStep)%io.esave==0)
                fwrite(ener,sizeof(double),1,efile);
//             fprintf(stderr,"\n");
            return MV_REJ ;
        }
    }
}

