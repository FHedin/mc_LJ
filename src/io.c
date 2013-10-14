#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "io.h"
#include "tools.h"

static uint32_t dcd_header_empty=1;

void write_xyz(ATOM at[], DATA *dat, uint64_t when, FILE *outf)
{
    recentre(at,dat);
    
    uint32_t i=0;
    fprintf(outf,"%d\n#step %d\n",dat->natom,when);
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
    
    for (i=0;i<nat;i++)
      fscanf(inpf,"%s %lf %lf %lf\n",at[i].sym,&(at[i].x),&(at[i].y),&(at[i].z));
    
}

void write_dcd(ATOM at[], DATA *dat, uint64_t when)
{
    recentre(at,dat);
    
    uint32_t i=0;
    uint32_t input_integer[2]= {0};

    if (dcd_header_empty)
    {
        char corp[4]="CORD";

        uint32_t  ICNTRL[20]= {0};
        ICNTRL[0]=ICNTRL[3]=dat->nsteps/io.trsave;
        ICNTRL[1]=ICNTRL[2]=1;
        ICNTRL[19]=37;	//charmm version

        uint32_t NTITLE=2;
        char TITLE[80]="";

        uint32_t NATOM=dat->natom;

        input_integer[0] = sizeof(corp) + sizeof(ICNTRL);
        input_integer[1] = sizeof(NTITLE);

        fwrite(&input_integer[0],sizeof(uint32_t),1,traj);
        {
            fwrite(corp,sizeof(char),4,traj);
            fwrite(ICNTRL,sizeof(uint32_t),20,traj);
        }
        fwrite(&input_integer[0],sizeof(uint32_t),1,traj);


        fwrite(&input_integer[0],sizeof(uint32_t),1,traj);
        {
            fwrite(&NTITLE,sizeof(uint32_t),1,traj);
            for (i=0; i<NTITLE; i++)
                fwrite(TITLE,sizeof(char),80,traj);
        }
        fwrite(&input_integer[1],sizeof(uint32_t),1,traj);


        fwrite(&input_integer[1],sizeof(uint32_t),1,traj);
        fwrite(&NATOM,sizeof(uint32_t),1,traj);
        fwrite(&input_integer[1],sizeof(uint32_t),1,traj);

        dcd_header_empty=0;
    }

    float x=0.f,y=0.f,z=0.f;
    input_integer[1]=sizeof(float)*dat->natom;

    fwrite(&input_integer[1],sizeof(uint32_t),1,traj);
    for(i=0; i<dat->natom; i++)
    {
        x=(float)at[i].x;
        fwrite(&x,sizeof(float),1,traj);
    }
    fwrite(&input_integer[1],sizeof(uint32_t),1,traj);

    fwrite(&input_integer[1],sizeof(uint32_t),1,traj);
    for(i=0; i<dat->natom; i++)
    {
        y=(float)at[i].y;
        fwrite(&y,sizeof(float),1,traj);
    }
    fwrite(&input_integer[1],sizeof(uint32_t),1,traj);

    fwrite(&input_integer[1],sizeof(uint32_t),1,traj);
    for(i=0; i<dat->natom; i++)
    {
        z=(float)at[i].z;
        fwrite(&z,sizeof(float),1,traj);
    }
    fwrite(&input_integer[1],sizeof(uint32_t),1,traj);

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