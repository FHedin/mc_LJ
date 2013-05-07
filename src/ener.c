#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "ener.h"

#ifndef K_CONSTRAINT
#define K_CONSTRAINT	2.50
#endif

/* How to call this function :
 *
 *  get_LJV(at,&dat,-1) is for total energy of the whole system.
 *
 *  get_LJV(at,&dat,candidate_atom_number) is for energy evaluation of candidate_atom_number only.
 *
 */
double get_LJ_V(ATOM at[],DATA *dat,int candidate)
{

    int i,j;
    double dx1,dy1,dz1;
    double dx2,dy2,dz2;
    double dcm;
    double d2, epsi_g, sig_g;
    double energy = 0.0;
    
    dat->E_constr = 0.0;
    CM cm = getCM(at,dat);

    if (candidate==-1)
    {
        for (i=0; i<(dat->natom); i++)
        {
            dx1=at[i].x;
            dy1=at[i].y;
            dz1=at[i].z;

            dcm = X2(cm.cx-dx1) +  X2(cm.cy-dy1) + X2(cm.cz-dz1) ;
            dat->E_constr += getExtraPot(dcm,at[i].ljp.sig,at[i].ljp.eps);

            for (j=i+1; j<(dat->natom); j++)
            {
                dx2=at[j].x;
                dy2=at[j].y;
                dz2=at[j].z;

                d2 = X2(dx2-dx1) +  X2(dy2-dy1) + X2(dz2-dz1) ;
                epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
                sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );

                energy += 4.0 * epsi_g *( (X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) );
//                 dat->E_constr += getExtraPot(d2,sig_g,epsi_g);
            }
        }
    }
    else
    {
        i=candidate;

        dx1=at[i].x;
        dy1=at[i].y;
        dz1=at[i].z;

        dcm = X2(cm.cx-dx1) +  X2(cm.cy-dy1) + X2(cm.cz-dz1) ;
        dat->E_constr += getExtraPot(dcm,at[i].ljp.sig,at[i].ljp.eps);

        for (j=0; j<(dat->natom); j++)
        {
            if (j!=i)
            {
                dx2=at[j].x;
                dy2=at[j].y;
                dz2=at[j].z;

                d2 = X2(dx2-dx1) +  X2(dy2-dy1) + X2(dz2-dz1) ;
                epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
                sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );

                energy += 4.0 * epsi_g *( (X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) );
//            	dat->E_constr += getExtraPot(d2,sig_g,epsi_g);
            }
        }
    }

    return energy;
}

void get_LJ_DV(ATOM at[], DATA *dat, double DV[])
{
    //DV[] is 3*natom length : first x values, then y then z.

    int i=0 , j=0 ;
    double epsi_g=0.0 , sig_g=0.0;
    double dx=0.0 , dy=0.0 , dz=0.0 , d2=0.0 ;
    double DE=0.0 ;

    for (i=0 ; i < dat->natom ; i++ )
    {
        DV[i] = 0.0 ;
        DV[i+dat->natom] = 0.0 ;
        DV[i+2*dat->natom] = 0.0 ;
        for (j=0 ; j < dat->natom ; j++ )
        {
            if (i==j) continue ;
            epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
            sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );
            dx = at[i].x - at[j].x ;
            dy = at[i].y - at[j].y ;
            dz = at[i].z - at[j].z ;
            d2  = dx*dx + dy*dy + dz*dz ;
            DE = -24.0*epsi_g*( 2.0*(X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) ) / d2 ;
            DV[i] += DE*dx;
            DV[i+dat->natom] += DE*dy;
            DV[i+2*dat->natom] += DE*dz;
        }
    }
}

double get_AZIZ_V(ATOM at[], DATA *dat, int candidate)
{

    int i,j;
    double d, d2, energy=0.0, etmp;

    if (candidate==-1)
    {
        for (i=0; i<(dat->natom-1); i++)
        {
            for (j=i+1; j<(dat->natom); j++)
            {
                d2 = X2(at[i].x-at[j].x) +  X2(at[i].y-at[j].y) + X2(at[i].z-at[j].z) ;
                d = sqrt(d2);

                if(!strcmp(at[i].sym,at[j].sym)) //if the same type (strcmp return 0 if identical)
                {
                    if (!strcmp(at[i].sym,"Ne")) //both are neon
                        aziz_ne_ne_(&d,&etmp);
                    else	//both are argon
                        aziz_ar_ar_(&d,&etmp);

                }
                else //if different type it means it is 1 Ar and 1 Ne
                    aziz_ar_ne_(&d,&etmp);

                energy += etmp;
            }
        }
    }
    else
    {
        i=candidate;
        for (j=0; j<(dat->natom); j++)
        {
            if (j!=i)
            {
                d2 = X2(at[i].x-at[j].x) +  X2(at[i].y-at[j].y) + X2(at[i].z-at[j].z) ;
                d = sqrt(d2);

                if(!strcmp(at[i].sym,at[j].sym)) //if the same type (strcmp return 0 if identical)
                {
                    if (!strcmp(at[i].sym,"Ne")) //both are neon
                        aziz_ne_ne_(&d,&etmp);
                    else	//both are argon
                        aziz_ar_ar_(&d,&etmp);

                }
                else //if different type it means it is 1 Ar and 1 Ne
                    aziz_ar_ne_(&d,&etmp);

                energy += etmp;
            }
        }
    }

    return energy*CM1TOKJM*JTOCAL;
}

double getExtraPot(double d2, double sig, double eps)
{
	double vc = d2/(X2(K_CONSTRAINT*sig));
	vc=pow(vc,10.0);
	vc*=eps;

	return vc;
}

