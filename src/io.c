/**
 * \file io.c
 *
 * \brief Input/Output related functions
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
#include <time.h>
#include <string.h>

#include "global.h"
#include "io.h"
#include "tools.h"

/// header has to be written only once at the beginning of the dcd
static uint32_t dcd_header_empty=1;

/**
 * Writes one xyz file
 * 
 * @param at ATOM array where to store coordinates
 * @param dat common simulation data
 * @param when Step at which we save data
 * @param outf FILE where to save xyz
 */
void write_xyz(ATOM at[], DATA *dat, uint64_t when, FILE *outf)
{
    recentre(at,dat);

    uint32_t i=0;
    fprintf(outf,"%d\n#step %"PRIu64"\n",dat->natom,when);
    for (i=0; i<(dat->natom); i++)
        fprintf(outf,"%s\t%10.5lf\t%10.5lf\t%10.5lf\n",at[i].sym,at[i].x,at[i].y,at[i].z);
}

/**
 * Reads one xyz file
 * 
 * @param at ATOM array where to store coordinates
 * @param dat common simulation data
 * @param inpf FILE where to read xyz from
 */
void read_xyz(ATOM at[], DATA *dat, FILE *inpf)
{
    uint32_t i;
    uint32_t nat=0;
    char comm[1024]="";

    fscanf(inpf,"%d\n",&nat);
    fgets(comm,1024,inpf);

    if (nat != dat->natom)
    {
        fprintf(stdout,"Error : natom in the xyz file differs from the declaration in the input file.\n");
        exit(-1);
    }

    for (i=0; i<nat; i++)
        fscanf(inpf,"%s %lf %lf %lf\n",at[i].sym,&(at[i].x),&(at[i].y),&(at[i].z));

}

/**
 * Writes a CHARMM like dcd
 * @param at ATOM array where to store coordinates
 * @param dat common simulation data
 * @param when At which step function was called
 */
void write_dcd(ATOM at[], DATA *dat, uint64_t when)
{
    recentre(at,dat);

    uint32_t i=0;
    uint32_t sizeB = 0;

    if (dcd_header_empty)
    {
        char corp[4]= {'C','O','R','D'};

        uint32_t  ICNTRL[20]= {0};
        ICNTRL[0]=ICNTRL[3] = (uint32_t) dat->nsteps/io.trsave;
        ICNTRL[1]=ICNTRL[2]=1;
        ICNTRL[19]=39;	//charmm version

        uint32_t NTITLE=3;
        char TITLE[3][80]= {"","",""};

        // add to title info about date, user, machine, etc... if available
        time_t rawtime;
        struct tm *timeinfo;
        char *user=NULL , *host=NULL , *pwd=NULL;
        char tmp[4096]="";
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        user=getenv("USER");
        host=getenv("HOSTNAME");
        pwd=getenv("PWD");
        sprintf(tmp,"* CREATION DATE : %s",asctime(timeinfo));
        // strncat : char * strncat ( char * destination, const char * source, size_t num );
        // Appends the first num characters of source to destination, plus a terminating null-character.
        strncat(TITLE[0],tmp,79);
        sprintf(tmp,"* USER : %s HOSTNAME : %s",(user!=NULL)?user:"UNKNOWN",(host!=NULL)?host:"UNKNOWN");
        strncat(TITLE[1],tmp,79);
        sprintf(tmp,"* PWD : %s",(pwd!=NULL)?pwd:"UNKNOWN");
        strncat(TITLE[2],tmp,79);

        uint32_t NATOM=dat->natom;

        sizeB = sizeof(corp) + sizeof(ICNTRL);
        fwrite(&sizeB,sizeof(uint32_t),1,traj);
        {
            fwrite(corp,sizeof(char),4,traj);
            fwrite(ICNTRL,sizeof(uint32_t),20,traj);
        }
        fwrite(&sizeB,sizeof(uint32_t),1,traj);

        sizeB = sizeof(NTITLE) + NTITLE*80*sizeof(char);
        fwrite(&sizeB,sizeof(uint32_t),1,traj);
        {
            fwrite(&NTITLE,sizeof(uint32_t),1,traj);
            for (i=0; i<NTITLE; i++)
                fwrite(TITLE[i],sizeof(char),80,traj);
        }
        fwrite(&sizeB,sizeof(uint32_t),1,traj);

        sizeB = sizeof(NATOM);
        fwrite(&sizeB,sizeof(uint32_t),1,traj);
        fwrite(&NATOM,sizeof(uint32_t),1,traj);
        fwrite(&sizeB,sizeof(uint32_t),1,traj);

        dcd_header_empty=0;
    }

    float x=0.f,y=0.f,z=0.f;
    sizeB=(uint32_t)sizeof(float)*dat->natom;

    fwrite(&sizeB,sizeof(uint32_t),1,traj);
    for(i=0; i<dat->natom; i++)
    {
        x=(float)at[i].x;
        fwrite(&x,sizeof(float),1,traj);
    }
    fwrite(&sizeB,sizeof(uint32_t),1,traj);

    fwrite(&sizeB,sizeof(uint32_t),1,traj);
    for(i=0; i<dat->natom; i++)
    {
        y=(float)at[i].y;
        fwrite(&y,sizeof(float),1,traj);
    }
    fwrite(&sizeB,sizeof(uint32_t),1,traj);

    fwrite(&sizeB,sizeof(uint32_t),1,traj);
    for(i=0; i<dat->natom; i++)
    {
        z=(float)at[i].z;
        fwrite(&z,sizeof(float),1,traj);
    }
    fwrite(&sizeB,sizeof(uint32_t),1,traj);

}
/**
 * Writes a restart file : DO NOT USE for the moment there is a bug somewhere !
 * 
 * @bug This function is bugged for the moment please don't use it !!
 * 
 * @param at Atom array
 * @param dat Data structure
 * @param spdat Spatial Averaging Data structure
 * @param meth method used : 0 = metropolis and 1 = SA-MC
 */
void write_rst(ATOM at[], DATA *dat, SPDAT *spdat, uint32_t meth)
{
    FILE *rstfile=NULL;
    rstfile=fopen("restart.dat","wb");

    /** From global.h **/
    // 1 : write global variables
    fwrite(&charmm_units,sizeof(uint32_t),1,rstfile);

    // 2 : structure DATA
    fwrite(&dat->natom,sizeof(uint32_t),1,rstfile);
    fwrite(dat->method,sizeof(char),32,rstfile);
    fwrite(&dat->nsteps,sizeof(uint64_t),1,rstfile);

    fwrite(&dat->d_max,sizeof(double),1,rstfile);
    fwrite(&dat->d_max_when,sizeof(uint32_t),1,rstfile);
    fwrite(&dat->d_max_tgt,sizeof(double),1,rstfile);

    fwrite(&dat->inid,sizeof(double),1,rstfile);
    fwrite(&dat->T,sizeof(double),1,rstfile);
    fwrite(&dat->E_steepD,sizeof(double),1,rstfile);
    fwrite(&dat->E_expected,sizeof(double),1,rstfile);
    fwrite(&dat->beta,sizeof(double),1,rstfile);

    if(meth==1)
    {
        // 2 : structure SPDAT
        fwrite(&spdat->meps,sizeof(uint32_t),1,rstfile);
        fwrite(&spdat->neps,sizeof(uint32_t),1,rstfile);
        fwrite(&spdat->weps,sizeof(double),1,rstfile);
    }

    // 3 : structure ATOM
    for(uint32_t i=0; i<dat->natom; i++)
    {
        fwrite(&(at[i].x),sizeof(double),1,rstfile);
        fwrite(&(at[i].y),sizeof(double),1,rstfile);
        fwrite(&(at[i].z),sizeof(double),1,rstfile);
        
        fwrite(at[i].sym,sizeof(char),4,rstfile);
        
        fwrite(at[i].ljp.sym,sizeof(char),4,rstfile);
        
        fwrite(&(at[i].ljp.sig),sizeof(double),1,rstfile);
        fwrite(&(at[i].ljp.eps),sizeof(double),1,rstfile);
    }

    /** END **/
    fclose(rstfile);
}
