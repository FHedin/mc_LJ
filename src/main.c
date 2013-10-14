#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if __unix__
#include <sys/resource.h>
#endif /* for some unixes */

#include "global.h"
#include "MCclassic.h"
#include "MCspav.h"
#include "tools.h"
#include "rand.h"
#include "ener.h"
#include "io.h"
#include "parsing.h"

//define where is the null file
#if __unix__
#define NULLFILE "/dev/null"
#elif _WIN32
#define NULLFILE "nul"
#endif

/**
 * 	Some global variables
 **/

/*
 *     char crdtitle_first[128],char crdtit*le_last[128],char trajtitle[128],char etitle[128];
 *     uint32_t esave,uint32_t trsave;
 */
IODAT io = {NULLFILE,NULLFILE,NULLFILE,NULLFILE,1000,1000};

//common use files
FILE *traj=NULL;
FILE *crdfile=NULL;
FILE *efile=NULL;

//is the stdout redirected ?
uint32_t is_stdout_redirected=0;
//are we using charmm units ?
uint32_t charmm_units=0;

#ifdef _OPENMP
//for para execution
uint32_t ncpus=1,nthreads=1;
#endif

/**
 * End
 **/

//prototypes
int  main(int argc, char **argv);
void start_classic(DATA *dat, ATOM at[]);
void start_spav(DATA *dat, SPDAT *spdat, ATOM at[]);
void help(char **argv);
void getValuesFromDB(DATA *dat);

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        fprintf(stdout,"Error : no input file.\n");
        help(argv);
        return EXIT_SUCCESS;
    }

    //os-independant redirection of stderr to the null file
    freopen(NULLFILE,"w",stderr);

    uint32_t i;
    char seed[128] = "";
    char inpf[128] = "";

    DATA dat ;
    SPDAT spdat = {5,5,0.5,NULL,0};
    ATOM *at = NULL;

    get_ENER = NULL;
    get_DV = NULL;
    write_traj= &(write_dcd);

    //arguments parsing
    for (i=1; i<(uint32_t)argc; i++)
    {
        if (!strcasecmp(argv[i],"-i"))
            sprintf(inpf,"%s",argv[++i]);
        else if (!strcasecmp(argv[i],"-seed"))
            sprintf(seed,"%s",argv[++i]);
        else if (!strcasecmp(argv[i],"-of"))
        {
            freopen(argv[++i],"w",stdout);
            is_stdout_redirected = 1 ;
        }
        else if (!strcasecmp(argv[i],"-ef"))
        {
            freopen(argv[++i],"w",stderr);
        }
#ifdef _OPENMP
        else if (!strcasecmp(argv[i],"-np"))
        {
            nthreads=atoi(argv[++i]);
            ncpus=omp_get_num_procs();
            (nthreads < ncpus) ? omp_set_num_threads(nthreads) : omp_set_num_threads(ncpus);
            omp_set_num_threads(nthreads);
        }
#endif
        else
        {
            fprintf(stdout,"Error : Argument '%s' not recognised.\n",argv[i]);
            exit(-2);
        }
    }

#ifdef STDRAND
    srand(time(NULL));
    dat.nrn = 2048 ;
    dat.rn = calloc(dat.nrn,sizeof dat.rn);
#else
    if (!strlen(seed))
        sprintf(seed,"%d",(int)time(NULL)) ;
    dat.nrn = 2048 ;
    dat.rn = calloc(dat.nrn,sizeof dat.rn);
    dat.seeds = calloc(strlen(seed),sizeof dat.seeds);
    for (i=0; i<(uint32_t)strlen(seed); i++)
    {
        dat.seeds[i]=(uint32_t)seed[strlen(seed)-1-i];
        dat.seeds[i]*=dat.seeds[0];
    }
    dsfmt_init_by_array(&(dat.dsfmt),dat.seeds,(int32_t)strlen(seed));
#endif

    at=parse_from_file(inpf,&dat,&spdat);

    //set the pointer to the default energy function
    if(get_ENER==NULL && get_DV==NULL)
    {
        get_ENER = &(get_LJ_V);
        get_DV = &(get_LJ_DV);
    }

#ifdef _OPENMP
    fprintf(stdout,"Starting program with %d threads (%d cpus available)\n\n",nthreads,ncpus);
#else
    fprintf(stdout,"Starting program in sequential mode\n\n");
#endif

    fprintf(stdout,"SEED   = %s \n\n",seed);

    if (get_ENER==&(get_LJ_V))          fprintf(stdout,"Using L-J potential\n");
    else if (get_ENER==&(get_AZIZ_V))   fprintf(stdout,"Using Aziz potential\n");

    if (charmm_units)                   fprintf(stdout,"Using CHARMM  units.\n\n");
    else				fprintf(stdout,"Using REDUCED units.\n\n");

    fprintf(stdout,"Energy      saved each %d  steps in file %s\n",io.esave,io.etitle);
    fprintf(stdout,"Trajectory  saved each %d  steps in file %s\n",io.trsave,io.trajtitle);
    fprintf(stdout,"Initial configuration saved in file %s\n",io.crdtitle_first);
    fprintf(stdout,"Final   configuration saved in file %s\n\n",io.crdtitle_last);
    
    getValuesFromDB(&dat);

    if (charmm_units)
        dat.beta = 1.0/(KBCH*dat.T);
    else
        dat.beta = 1.0/(dat.T);

    fprintf(stdout,"method = %s\n",dat.method);
    fprintf(stdout,"natom  = %d\n",dat.natom);
    fprintf(stdout,"nsteps = %ld\n",dat.nsteps);
    fprintf(stdout,"T      = %lf \n",dat.T);
    fprintf(stdout,"beta   = %lf\n",dat.beta);
    
    if(dat.d_max_when==0)
        fprintf(stdout,"dmax   = %lf (fixed) \n\n",dat.d_max);
    else
        fprintf(stdout,"dmax   = %4.2lf updated each %d steps for "
                       "targeting %4.2lf %% of acceptance \n\n",dat.d_max,dat.d_max_when,dat.d_max_tgt);

    if (strcasecmp(dat.method,"metrop")==0)
    {
        start_classic(&dat,at);
    }
    else if (strcasecmp(dat.method,"spav")==0)
    {
        start_spav(&dat,&spdat,at);
    }
    else
    {
        fprintf(stderr,"Error : Method [%s] unknowm.\n",dat.method);
        free(dat.rn);
#ifndef STDRAND
        free(dat.seeds);
#endif
        exit(-3);
    }

#if ( defined __linux__ || defined __FreeBSD__ )
    //compatible with some unixes-like OS: the struct rusage communicates with the kernel directly.
    struct rusage infos_usage;
    getrusage(RUSAGE_SELF,&infos_usage);
    fprintf(stdout,"Memory used in kBytes is : %ld\n",infos_usage.ru_maxrss);
    fprintf(stdout,"Execution time in Seconds : %lf\n",
           (double)infos_usage.ru_utime.tv_sec+(double)infos_usage.ru_utime.tv_usec/1000000.0);
#endif /* unixes part */
    fprintf(stdout,"End of program\n");
    
    //write a rst file
    //write_rst(&dat,&spdat,at);

    free(dat.rn);
#ifndef STDRAND
    free(dat.seeds);
#endif

    free(at);

    return EXIT_SUCCESS;
}

void start_classic(DATA *dat, ATOM at[])
{
    double ener = 0.0 ;
    uint64_t acc=0;

    crdfile=fopen(io.crdtitle_first,"wt");
    efile=fopen(io.etitle,"wb");
    traj=fopen(io.trajtitle,"wb");
    
    write_xyz(at,dat,0,crdfile);
    fclose(crdfile);
    
    // 	sprintf(trname,"%d.dcd",i+1);

    freopen(io.trajtitle,"wb",traj);

    ener = (*get_ENER)(at,dat,-1);

    fprintf(stdout,"\nStarting METROP Monte-Carlo\n");
    fprintf(stdout,"LJ initial energy is : %lf \n\n",ener);
    fwrite(&ener,sizeof(double),1,efile);

    acc=make_MC_moves(at,dat,&ener);

    fprintf(stdout,"\n\nLJ final energy is : %lf\n",ener);
    fprintf(stdout,"Acceptance ratio is %lf %% \n",100.0*(double)acc/(double)dat->nsteps);
    fprintf(stdout,"Final dmax = %lf\n",dat->d_max);
    fprintf(stdout,"End of METROP Monte-Carlo\n\n");

    crdfile=fopen(io.crdtitle_last,"wt");
    write_xyz(at,dat,dat->nsteps,crdfile);
    
    fclose(crdfile);
    fclose(traj);
    fclose(efile);

}

void start_spav(DATA *dat, SPDAT *spdat, ATOM at[])
{
    fprintf(stdout,"SPAV parameters are :\n");
    fprintf(stdout,"W_EPSILON = %lf\nM_EPSILON = %d\nN_EPSILON = %d\n\n",spdat->weps,spdat->meps,spdat->neps);

    spdat->normalSize=2048;
    spdat->normalNumbs=malloc(spdat->normalSize*sizeof spdat->normalNumbs);

    double ener = 0.0 ;
    uint64_t acc=0;

    crdfile=fopen(io.crdtitle_first,"wt");
    efile=fopen(io.etitle,"wb");
    traj=fopen(io.trajtitle,"wb");
    
    write_xyz(at,dat,0,crdfile);
    fclose(crdfile);
    
    ener = (*get_ENER)(at,dat,-1);

    fprintf(stdout,"\nStarting SPAV\n");
    fprintf(stdout,"LJ initial energy is : %lf \n\n",ener);
    fwrite(&ener,sizeof(double),1,efile);

    acc=launch_SPAV(at,dat,spdat,&ener);

    fprintf(stdout,"LJ final energy is : %lf\n",ener);
    fprintf(stdout,"Acceptance ratio is %lf %% \n",100.0*(double)acc/(double)dat->nsteps);
    fprintf(stdout,"final dmax = %lf\n",dat->d_max);
    fprintf(stdout,"End of SPAV\n\n");
    
    crdfile=fopen(io.crdtitle_last,"wt");
    write_xyz(at,dat,dat->nsteps,crdfile);

    fclose(crdfile);
    fclose(traj);
    fclose(efile);

    // 	fclose(stfile);

    free(spdat->normalNumbs);
}

void help(char **argv)
{
    fprintf(stdout,"Need an argument : %s -i an_input_file\n",argv[0]);
    fprintf(stdout,"{optional args : -seed [a_rnd_seed] -of [output_file] -ef [error_file]}\nExample:\n");
    fprintf(stdout,"%s -i input_file -seed 1330445520 -of out.log -ef err.log \n",argv[0]);
}


void getValuesFromDB(DATA *dat)
{
    if (dat->natom==13)
    {
        dat->E_steepD   = -33.0;
        dat->E_expected = -44.326801;
    }
    else if (dat->natom==19)
    {
        dat->E_steepD   = -66.0;
        dat->E_expected = -72.659782;
    }
    else if (dat->natom==31)
    {
        dat->E_steepD   = -129.0;
        dat->E_expected = -133.586422;
    }
    else if (dat->natom==37)
    {
        dat->E_steepD   = -162;
        dat->E_expected = -167.033672;
    }
    else if (dat->natom==38)
    {
        dat->E_steepD   = -167.5;
        dat->E_expected = -173.928427;
    }
    else if (dat->natom==55)
    {
        dat->E_steepD   = -270.5;
        dat->E_expected = -279.248470;
    }
    else if (dat->natom==75)
    {
        dat->E_steepD   = -388.0;
        dat->E_expected = -397.492331;
    }
    else
    {
        dat->E_steepD   = -999999999.999999;
        dat->E_expected = -999999999.999999;
    }
}

