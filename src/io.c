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
#include <time.h>
#include <string.h>

#include "global.h"
#include "io.h"
#include "tools.h"

static uint32_t dcd_header_empty=1;

void write_xyz(ATOM at[], DATA *dat, uint64_t when, FILE *outf)
{
    recentre(at,dat);

    uint32_t i=0;
    fprintf(outf,"%d\n#step %"PRIu64"\n",dat->natom,when);
    for (i=0; i<(dat->natom); i++)
        fprintf(outf,"%s\t%10.5lf\t%10.5lf\t%10.5lf\n",at[i].sym,at[i].x,at[i].y,at[i].z);
}

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

void write_rst(DATA *dat, SPDAT *spdat, ATOM at[])
{
    FILE *rstfile=NULL;
    rstfile=fopen("restart.dat","wb");

    /** From global.h **/
    // 1 : write global variables
    fwrite(&is_stdout_redirected,sizeof(uint32_t),1,rstfile);
    fwrite(&charmm_units,sizeof(uint32_t),1,rstfile);
#ifdef _OPENMP
    fwrite(&ncpus,sizeof(uint32_t),1,rstfile);
    fwrite(&nthreads,sizeof(uint32_t),1,rstfile);
#endif

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

    // 2 : structure SPDAT
    fwrite(&spdat->meps,sizeof(uint32_t),1,rstfile);
    fwrite(&spdat->neps,sizeof(uint32_t),1,rstfile);
    fwrite(&spdat->weps,sizeof(double),1,rstfile);
    fwrite(&spdat->normalSize,sizeof(uint32_t),1,rstfile);
    fwrite(spdat->normalNumbs,sizeof(double),spdat->normalSize,rstfile);

    // 3 : structure ATOM
    fwrite(at,sizeof(ATOM),dat->natom,rstfile);

    /** END **/
    //fclose(rstfile);
}
