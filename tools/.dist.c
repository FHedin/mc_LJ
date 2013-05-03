#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define X2(x) (x)*(x)

/**
**	This programs estimates the distance between the first aArgon atom and the center of mass of all the other Neon atoms
**/

typedef struct
{
	double x,y,z ;
	char sym[4] ;
}ATOM;

typedef struct
{
	double x,y,z ;
}CENTER;

void read_onestep_xyz(FILE *tr, ATOM *at, int natom);

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		puts("Error : trajectory file expected");
		return -1;
	}
	
	FILE *tr = NULL ;
	tr = fopen(argv[1],"r") ;
	
	FILE *out = NULL;
	out = fopen("d_ar.dat","w");
	
	ATOM *at = NULL ;
	CENTER ct = {0.,0.,0.};
	
	int natom=0;
	char comm_line[32]="";
	int st=0;
	
	int i=0;
	double d=0.;
	
	fscanf(tr,"%d",&natom);
	at=malloc(natom*sizeof(ATOM));
	
	do
	{
		fscanf(tr,"%s %d",comm_line,&st);
		read_onestep_xyz(tr,at,natom);
		for (i=1;i<natom;i++)
		{
			ct.x += at[i].x;
			ct.y += at[i].y;
			ct.z += at[i].z;
		}
		ct.x /= natom-1;
		ct.y /= natom-1;
		ct.z /= natom-1;
		d = sqrt ( X2(at[0].x-ct.x) + X2(at[0].y-ct.y) + X2(at[0].z-ct.z) ) ;
		fprintf(out,"%lf\n",d);
		if(fscanf(tr,"%d",&natom)==EOF)
			break;
	}while(1);
	
}

void read_onestep_xyz(FILE *tr, ATOM *at, int natom)
{
	int i=0;
	for (i=0;i<natom;i++)
		fscanf(tr,"%s %lf %lf %lf",at[i].sym,&(at[i].x),&(at[i].y),&(at[i].z));
}

