#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "io.h"

void write_xyz(ATOM at[], DATA *dat, int when)
{
	int i=0;
	fprintf(traj,"%d\n#step %d\n",dat->natom,when);
	for (i=0; i<(dat->natom); i++)
		fprintf(traj,"%s\t%10.5lf\t%10.5lf\t%10.5lf\n",at[i].sym,at[i].x,at[i].y,at[i].z);
}

void write_dcd(ATOM at[], DATA *dat, int when)
{
	unsigned int i=0;
	unsigned int input_integer[2]={0};
	
	if (dcd_header_empty)
	{
		char corp[4]="CORD";
		
		unsigned int  ICNTRL[20]={0};
		ICNTRL[0]=ICNTRL[3]=((dat->nsteps/io.trsave)+1);
		ICNTRL[1]=ICNTRL[2]=1;
		ICNTRL[19]=37;	//charmm version
		
		unsigned int NTITLE=2;
		char TITLE[80]="";
		
		unsigned int NATOM=dat->natom;
		
		input_integer[0] = sizeof(corp) + sizeof(ICNTRL);
		input_integer[1] = sizeof(NTITLE);
		
		fwrite(&input_integer[0],sizeof(unsigned int),1,traj);
		{
			fwrite(corp,sizeof(char),4,traj);
			fwrite(ICNTRL,sizeof(int),20,traj);
		}
		fwrite(&input_integer[0],sizeof(unsigned int),1,traj);
		

		fwrite(&input_integer[0],sizeof(unsigned int),1,traj);
		{
			fwrite(&NTITLE,sizeof(int),1,traj);
			for (i=0;i<NTITLE;i++)
				fwrite(TITLE,sizeof(char),80,traj);
		}
		fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
		
		
		fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
			fwrite(&NATOM,sizeof(int),1,traj);
		fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
		
		dcd_header_empty=0;
	}
	
	float x=0.f,y=0.f,z=0.f;
	input_integer[1]=sizeof(float)*dat->natom;
	
	fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
	for(i=0;i<dat->natom;i++)
	{
		x=(float)at[i].x;
		fwrite(&x,sizeof(float),1,traj);
	}
	fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
	
	fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
	for(i=0;i<dat->natom;i++)
	{
		y=(float)at[i].y;
		fwrite(&y,sizeof(float),1,traj);
	}
	fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
	
	fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
	for(i=0;i<dat->natom;i++)
	{
		z=(float)at[i].z;
		fwrite(&z,sizeof(float),1,traj);
	}
	fwrite(&input_integer[1],sizeof(unsigned int),1,traj);
	
}

