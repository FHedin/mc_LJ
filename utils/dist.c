#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define X2(x) (x)*(x)

/**
**	This program estimates the distance between the first Argon atom and the center of mass of all the other Neon atoms
**/

typedef struct
{
    float x,y,z ;
} ATOM;

typedef struct
{
    float x,y,z ;
} CENTER;

void read_dcd_header(FILE *dcdf, int *natom, int *steps);
void read_onestep_dcd(FILE *dcdf, ATOM *at, int natom);

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        puts("Error : trajectory file expected");
        return EXIT_FAILURE;
    }

    FILE *tr = NULL ;
    tr = fopen(argv[1],"rb") ;

    FILE *out = NULL;
    out = fopen("d_ar.dat","w");

    ATOM *at = NULL ;
    CENTER ct = {0.f,0.f,0.f};

    int natom=0,steps=0;
    read_dcd_header(tr,&natom,&steps);
    at=malloc(natom*sizeof(ATOM));

    int   i=0,st=0;
    float d=0.f;

    for (st=0; st<steps; st++)
    {
        read_onestep_dcd(tr,at,natom);
        for (i=1; i<natom; i++)
        {
            ct.x += at[i].x;
            ct.y += at[i].y;
            ct.z += at[i].z;
        }
        ct.x /= natom-1;
        ct.y /= natom-1;
        ct.z /= natom-1;
        d = sqrt ( X2(at[0].x-ct.x) + X2(at[0].y-ct.y) + X2(at[0].z-ct.z) ) ;
        fprintf(out,"%f\n",d);
        printf("Step %10d done\r",st);
        fflush(stdout);
    }

    free(at);
    fclose(tr);
    fclose(out);

    printf("\n");

    return EXIT_SUCCESS;
}

void read_dcd_header(FILE *dcdf, int *natom, int *steps)
{
    int i=0;

    char HDR[5]="";
    int  ICNTRL[20]= {0};
    fseek(dcdf,4,SEEK_CUR);
    {
        fread(HDR,sizeof(char),4,dcdf);
        fread(ICNTRL,sizeof(int),20,dcdf);
    }
    fseek(dcdf,4,SEEK_CUR);

    int NTITLE=0;
    char TITLE[80]="";
    fseek(dcdf,4,SEEK_CUR);
    {
        fread(&NTITLE,sizeof(int),1,dcdf);
        for (i=0; i<NTITLE; i++)
            fread(TITLE,sizeof(char),80,dcdf);
    }
    fseek(dcdf,4,SEEK_CUR);

    int NATOM=0;
    fseek(dcdf,4,SEEK_CUR);
    {
        fread(&NATOM,sizeof(int),1,dcdf);
    }
    fseek(dcdf,4,SEEK_CUR);

    *natom = NATOM;
    *steps = ICNTRL[0];
}

void read_onestep_dcd(FILE *dcdf, ATOM *at, int natom)
{
    static int i=0;
    for(i=0; i<natom; i++)
    {
        fseek(dcdf,4,SEEK_CUR);
        {
            for(i=0; i<natom; i++)
                fread(&at[i].x,sizeof(float),1,dcdf);
        }
        fseek(dcdf,4,SEEK_CUR);

        fseek(dcdf,4,SEEK_CUR);
        {
            for(i=0; i<natom; i++)
                fread(&at[i].y,sizeof(float),1,dcdf);
        }
        fseek(dcdf,4,SEEK_CUR);

        fseek(dcdf,4,SEEK_CUR);
        {
            for(i=0; i<natom; i++)
                fread(&at[i].z,sizeof(float),1,dcdf);
        }
        fseek(dcdf,4,SEEK_CUR);
    }
}

