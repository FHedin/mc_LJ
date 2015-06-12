/**
 * \file main.c
 *
 * \brief C program for MC simulations of Lennard Jones clusters.
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
#include <strings.h>
#include <time.h>
#include <math.h>

// if parallel execution is implemented by using openMP
#ifdef _OPENMP
#include <omp.h>
#endif

// for some unixes : usage info (total time, memory ...)
#ifdef __unix__
#include <sys/resource.h>
#endif

#include "global.h"
#include "MCclassic.h"
#include "MCspav.h"
#include "tools.h"
#include "rand.h"
#include "ener.h"
#include "minim.h"
#include "io.h"
#include "parsing.h"
#include "logger.h"
#include "plugins_lua.h"

// -----------------------------------------------------------------------------------------

/*
 * Some global variables initialisation here
 */

/*
 * Initialise the io structure containing file names,
 * by default everything is discarded to NULLFILE
 */
IODAT io = {NULLFILE,NULLFILE,NULLFILE,NULLFILE,1000,1000};

// pointer to FILE for trajectory, coordinates, energy
FILE *traj=NULL;
FILE *crdfile=NULL;
FILE *efile=NULL;

/*
 * boolean like values
 * is the stdout redirected ?
 * are we using charmm units ?
 */
uint32_t is_stdout_redirected=0;
uint32_t charmm_units=0;

#ifdef _OPENMP
//for para execution we will try to get the number of cpus and threads available
uint32_t ncpus=1,nthreads=1;
#endif

/*
 * Errors, warning, etc ... --> logging.
 * Default Level is LOG_WARNING, which means that everything which is at least
 * a warning is printed (it includes Errors also).
 * See logger.h for other possibilities.
 */
LOG_LEVELS LOG_SEVERITY = LOG_WARNING;

/*
 * Type of Lua plugin, i.e. PAIR or FFI (see "plugins_lua.h")
 */
#ifdef LUA_PLUGINS
LUA_PLUGIN_TYPE lua_plugin_type = PAIR;
#endif

/*
 *  End of global variables initialisation
 */

// -----------------------------------------------------------------------------------------

//prototypes of functions written in this main.c
void start_classic(DATA *dat, ATOM at[]);
void start_spav(DATA *dat, SPDAT *spdat, ATOM at[]);
void help(char **argv);
void getValuesFromDB(DATA *dat);

// -----------------------------------------------------------------------------------------
/**
 * \brief   Main entry of the program
 * \details The 2 variables \b argc and \b argv are used for extracting
 *          command line parameters.
 * \param   argc Number of arguments, at least one as \b argv contains at least the program's name.
 * \param   argv[] Array of character strings containing the arguments.
 * \return  On exit returns EXIT_SUCCESS, EXIT_FAILURE otherwise.
 */
int main(int argc, char** argv)
{
    /* arguments parsing, we need at least "prog_name -i an_input_file"
     * prints some more instructions if needed
     */
    if (argc < 3)
    {
        fprintf(stdout,"[Info] No input file ! \n");
        help(argv);
        return EXIT_SUCCESS;
    }

    uint32_t i;
    char seed[128] = "";
    char inpf[FILENAME_MAX] = "";

    DATA dat ;
    SPDAT spdat = {5,5,0.5,NULL,0};
    ATOM *at = NULL;

    // function pointers for energy and gradient, and trajectory
    get_ENER = NULL;
    get_DV = NULL;
    write_traj= &(write_dcd);

    // arguments parsing
    for (i=1; i<(uint32_t)argc; i++)
    {
        // get name of input file
        if (!strcasecmp(argv[i],"-i"))
        {
            sprintf(inpf,"%s",argv[++i]);
        }
        // get user specified seed, 128 characters max, keep it as a string for the moment
        else if (!strcasecmp(argv[i],"-seed"))
        {
            sprintf(seed,"%s",argv[++i]);
        }
        // reopen stdout to user specified file
        else if (!strcasecmp(argv[i],"-o"))
        {
            freopen(argv[++i],"w",stdout);
            is_stdout_redirected = 1 ;
        }
        // specify the logging level
        else if (!strcasecmp(argv[i],"-log"))
        {
            if (!strcasecmp(argv[++i],"no"))
            {
                LOG_SEVERITY = LOG_NOTHING;
            }
            else if (!strcasecmp(argv[i],"err"))
            {
                LOG_SEVERITY = LOG_ERROR;
            }
            else if (!strcasecmp(argv[i],"warn"))
            {
                LOG_SEVERITY = LOG_WARNING;
            }
            else if (!strcasecmp(argv[i],"info"))
            {
                LOG_SEVERITY = LOG_INFO;
            }
            else if (!strcasecmp(argv[i],"dbg"))
            {
                LOG_SEVERITY = LOG_DEBUG;
            }
            else
                fprintf(stdout,"[Warning] Unknown log level '%s' : default value used.\n\n",argv[i]);
        }
#ifdef _OPENMP
        // if compiled with openMP the user can specify a maximum number of cpus to use
        // check if it is not higher than the real amount of cpus
        else if (!strcasecmp(argv[i],"-np"))
        {
            nthreads=atoi(argv[++i]);
            ncpus=omp_get_num_procs();
            (nthreads < ncpus) ? omp_set_num_threads(nthreads) : omp_set_num_threads(ncpus);
            //omp_set_num_threads(nthreads);
        }
#endif
        // print help and proper exit
        else if ( !strcasecmp(argv[i],"-h") || !strcasecmp(argv[i],"-help") || !strcasecmp(argv[i],"--help") )
        {
            help(argv);
            return EXIT_SUCCESS;
        }
        // error if unknown command line option
        else
        {
            fprintf(stdout,"[Error] Argument '%s' is unknown.\n",argv[i]);
            help(argv);
            exit(-2);
        }
    }

    //prepare log files if necessary
    init_logfiles();

    // Print date and some env. variables
    fprintf(stdout,"Welcome to %s ! Command line arguments succesfully parsed, now intialising parameters...\n\n",argv[0]);
    fprintf(stdout,"Logging level is : %s : see the documentation to see which .log files are generated, and what they contain.\n\n",get_loglevel_string());
    fprintf(stdout,"Now printing some local informations : \n");
    fprintf(stdout,"DATE : %s\n",get_time());
    fprintf(stdout,"HOSTNAME : %s\n",getenv("HOSTNAME"));
    fprintf(stdout,"USER : %s\n",getenv("USER"));
    fprintf(stdout,"PWD : %s\n",getenv("PWD"));

    /*
     * Random numbers can be generated by using the standard functions from the C library (no guarantee on the quality)
     * or by using the dSFMT generator (default, extremely stable, no redudancy )
     * if STDRAND is defined we use the C library.
     *
     * Random numbers are stored in a double array dat.rn, of size dat.nrn
     *
     * If no string seed was passed by command line we generate one by using the unix timestamp
     *  -For STDRAND, this is directly used for srand()
     *  -For dSFMT, an array of integers generated by using the string seed is sent to dsfmt_init_by_array
     *
     */
    if (!strlen(seed))
        sprintf(seed,"%d",(uint32_t)time(NULL)) ;
    LOG_PRINT(LOG_INFO,"seed = %s \n",seed);
    dat.nrn = 2048 ;
    dat.rn = calloc(dat.nrn,sizeof dat.rn);
#ifdef STDRAND
    srand( (uint32_t)labs(atol(seed)) );
#else
    dat.seeds = calloc(strlen(seed),sizeof dat.seeds);
    for (i=0; i<(uint32_t)strlen(seed); i++)
        dat.seeds[i] = (uint32_t) seed[i] << 8;
    for (i=0; i<(uint32_t)strlen(seed); i++)
        dat.seeds[i] *= (dat.seeds[strlen(seed)-1]+i+1);

    dsfmt_init_by_array(&(dat.dsfmt),dat.seeds,(int32_t)strlen(seed));

    for (i=0; i<(uint32_t)strlen(seed); i++)
    {
        LOG_PRINT(LOG_INFO,"dat.seeds[%d] = %d \n",strlen(seed)-1-i,dat.seeds[strlen(seed)-1-i]);
    }
#endif
    
    // parse input file, initialise atom list
    parse_from_file(inpf,&dat,&spdat,&at);

    // set the pointer to the default V and dV functions
    if(get_ENER==NULL)
        get_ENER = &(get_LJ_V);

    if(get_DV==NULL)
        get_DV = &(get_LJ_DV);

    // allocate arrays used by energy minimisation function
    alloc_minim(&dat);

    // sum up parameters to output file

#ifdef _OPENMP
    fprintf(stdout,"\nStarting program with %d threads (%d cpus available)\n\n",nthreads,ncpus);
#else
    fprintf(stdout,"\nStarting program in sequential mode\n\n");
#endif

    fprintf(stdout,"Seed   = %s \n\n",seed);

    if (get_ENER==&(get_LJ_V))
        fprintf(stdout,"Using L-J potential\n");
    else if (get_ENER==&(get_AZIZ_V))
        fprintf(stdout,"Using Aziz potential\n");
#ifdef LUA_PLUGINS
    else if (get_ENER==&(get_lua_V))
        fprintf(stdout,"Using plugin pair potential\n");
    else if (get_ENER==&(get_lua_V_ffi))
        fprintf(stdout,"Using plugin ffi potential\n");
#endif
    
    if (charmm_units)
        fprintf(stdout,"Using CHARMM  units.\n\n");
    else
        fprintf(stdout,"Using REDUCED units.\n\n");

    fprintf(stdout,"Energy      saved each %d  steps in file %s\n",io.esave,io.etitle);
    fprintf(stdout,"Trajectory  saved each %d  steps in file %s\n",io.trsave,io.trajtitle);
    fprintf(stdout,"Initial configuration saved in file %s\n",io.crdtitle_first);
    fprintf(stdout,"Final   configuration saved in file %s\n\n",io.crdtitle_last);

    // get values of best minima for several well known LJ clusters
    getValuesFromDB(&dat);

    // set the inverse temperature depending of the type of units used
    if (charmm_units)
        dat.beta = 1.0/(KBCH*dat.T);
    else
        dat.beta = 1.0/(dat.T);

    // again print parameters
    fprintf(stdout,"method = %s\n",dat.method);
    fprintf(stdout,"natom  = %d\n",dat.natom);
    fprintf(stdout,"nsteps = %"PRIu64"\n",dat.nsteps);
    fprintf(stdout,"T      = %lf \n",dat.T);
    fprintf(stdout,"beta   = %lf\n",dat.beta);

    if(dat.d_max_when==0)
        fprintf(stdout,"dmax   = %lf (fixed) \n\n",dat.d_max);
    else
        fprintf(stdout,"dmax   = %4.2lf updated each %d steps for "
                "targeting %4.2lf %% of acceptance \n\n",dat.d_max,dat.d_max_when,dat.d_max_tgt);

    // then depending of the type of simulation run calculation
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
        LOG_PRINT(LOG_ERROR,"Method [%s] unknowm.\n",dat.method);
        free(dat.rn);
#ifndef STDRAND
        free(dat.seeds);
#endif
        exit(-3);
    }

#ifdef __unix__
    // compatible with some unixes-like OS: the struct rusage communicates with the kernel directly.
    struct rusage infos_usage;
    getrusage(RUSAGE_SELF,&infos_usage);
    fprintf(stdout,"Maximum amount of memory used in kBytes is : %ld\n",infos_usage.ru_maxrss);
    fprintf(stdout,"Execution time in Seconds : %lf\n",
            (double)infos_usage.ru_utime.tv_sec+(double)infos_usage.ru_utime.tv_usec/1000000.0 +
            (double)infos_usage.ru_stime.tv_sec+(double)infos_usage.ru_stime.tv_usec/1000000.0
           );
#endif
    fprintf(stdout,"End of program\n");

    // free memory and exit properly
    free(dat.rn);
#ifndef STDRAND
    free(dat.seeds);
#endif
    free(at);
    dealloc_minim();

#ifdef LUA_PLUGINS
    end_lua();
#endif //LUA_PLUGINS

    // closing log files is the last thing to do as errors may occur at the end
    close_logfiles();

    return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------------------
/**
 * \brief   This function starts a standard Metropolis Monte Carlo simulation.
 *
 * \details This function is first in charge of opening all the output (coordinates, trajectory and energy) files.\n
 *          Then the function \b #make_MC_moves starting the simulation is called.\n
 *          In the end it prints results, close the files and goes back to the function \b #main.
 *
 * \param   dat is a structure containing control parameters common to all simulations.
 * \param   at[] is an array of structures ATOM containing coordinates and other variables.
 */
void start_classic(DATA *dat, ATOM at[])
{
    double ener = 0.0 ;
    uint64_t acc=0;

    //open required files
    crdfile=fopen(io.crdtitle_first,"wt");
    efile=fopen(io.etitle,"wb");
    traj=fopen(io.trajtitle,"wb");

    //write initial coordinates
    write_xyz(at,dat,0,crdfile);
    fclose(crdfile);

    //get initial energy of whole system
    ener = (*get_ENER)(at,dat,-1);
    fprintf(stdout,"\nStarting METROP Monte-Carlo\n");
    fprintf(stdout,"LJ initial energy is : %lf \n\n",ener);

    //CALL TO MAIN mc FUNCTION
    acc=make_MC_moves(at,dat,&ener);
    //simulation finished here
    
    fprintf(stdout,"\n\nLJ final energy is : %lf\n",ener);
    fprintf(stdout,"Acceptance ratio is %lf %% \n",100.0*(double)acc/(double)dat->nsteps);
    fprintf(stdout,"Final dmax = %lf\n",dat->d_max);
    fprintf(stdout,"End of METROP Monte-Carlo\n\n");

    //write last coordinates
    crdfile=fopen(io.crdtitle_last,"wt");
    write_xyz(at,dat,dat->nsteps,crdfile);
    
    fclose(crdfile);
    fclose(traj);
    fclose(efile);
}

// -----------------------------------------------------------------------------------------
/**
 * \brief   This function starts a Spatial Averaging Monte Carlo  (SA-MC) simulation.
 *
 * \details This function is first in charge of opening all the output (coordinates, trajectory and energy) files.\n
 *          Then the function \b #launch_SPAV starting the simulation is called.\n
 *          In the end it prints results, close the files and goes back to the function \b #main.
 *
 * \param   dat is a structure containing control parameters common to all simulations.
 * \param   spdat is a structure containing control parameters dedicated to Spatial Averaging simulations.
 * \param   at[] is an array of structures ATOM containing coordinates and other variables.
 */
void start_spav(DATA *dat, SPDAT *spdat, ATOM at[])
{
    fprintf(stdout,"SPAV parameters are :\n");
    fprintf(stdout,"W_EPSILON = %lf\nM_EPSILON = %d\nN_EPSILON = %d\n\n",spdat->weps,spdat->meps,spdat->neps);

    // rand numbers related stuff
    spdat->normalSize=2048;
    spdat->normalNumbs=malloc(spdat->normalSize*sizeof spdat->normalNumbs);

    double ener = 0.0 ;
    uint64_t acc=0;
    
    alloc_SAMC(spdat);

    //open files
    crdfile=fopen(io.crdtitle_first,"wt");
    efile=fopen(io.etitle,"wb");
    traj=fopen(io.trajtitle,"wb");

    //write ini crds
    write_xyz(at,dat,0,crdfile);
    fclose(crdfile);

    //get E of whole system
    ener = (*get_ENER)(at,dat,-1);

    fprintf(stdout,"\nStarting SPAV\n");
    fprintf(stdout,"LJ initial energy is : %lf \n\n",ener);

    //run sp avg simulation
    acc=launch_SPAV(at,dat,spdat,&ener);

    fprintf(stdout,"LJ final energy is : %lf\n",ener);
    fprintf(stdout,"Acceptance ratio is %lf %% \n",100.0*(double)acc/(double)dat->nsteps);
    fprintf(stdout,"final dmax = %lf\n",dat->d_max);
    fprintf(stdout,"End of SPAV\n\n");

    //write final crds
    crdfile=fopen(io.crdtitle_last,"wt");
    write_xyz(at,dat,dat->nsteps,crdfile);

    dealloc_SAMC(spdat);

    fclose(crdfile);
    fclose(traj);
    fclose(efile);

    free(spdat->normalNumbs);
}

// -----------------------------------------------------------------------------------------
/**
 * \brief   This function simply prints a basic help message.
 *
 * \details If any of \b -h or \b -help or \b --help are provided on the command line this help message is printed.\n
 *          If no command line parameter is present this message is also printed.\n
 *          If an unknown command line parameter is present this message is also printed.
 *
 * \param   argv Simply the same array of command line parameters, from function \b #main.
 */
void help(char **argv)
{
    fprintf(stdout,"Need at least one argument : %s -i an_input_file\n",argv[0]);
    fprintf(stdout,"optional args : -seed [a_rnd_seed] -o [output_file] -log [logging level, one of { no | err | warn | info | dbg }] \n");
    fprintf(stdout,"Example : \n %s -i input_file -seed 1330445520 -o out.txt -log info \n\n",argv[0]);
    fprintf(stdout,"The default logging level is 'warn' \n");
}

// -----------------------------------------------------------------------------------------
/**
 * \brief   This sets variables possibly used when minimising the energy function.
 *
 * \details Depending of the value of \b dat->natom this sets the variables \b dat->E_steepD,
 *          which may be used as a threshold for starting Steepest Descent minimisation, and
 *          \b dat->E_expected which is the energy of the best minimum.
 *
 * \param   dat is a structure containing control parameters common to all simulations.
 */
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

// -----------------------------------------------------------------------------------------

