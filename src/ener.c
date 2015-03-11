/**
 * \file ener.c
 *
 * \brief Functions for getting hard coded Potential and Gradient (Lennard-Jones and Aziz)
 *          For used defined lua potentials see instead plugins.lua.c
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
#include <string.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "ener.h"

#ifndef K_CONSTRAINT
#define K_CONSTRAINT    4.00
#endif

/* How to call this function :
 *
 *  get_LJV(at,&dat,-1) is for total energy of the whole system.
 *
 *  get_LJV(at,&dat,candidate_atom_number) is for energy evaluation of candidate_atom_number only.
 *
 */
double get_LJ_V(ATOM at[], DATA *dat, int32_t candidate)
{

    uint32_t i,j;
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
            }
        }
    }
    else
    {
        i = (uint32_t) candidate;

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
            }
        }
    }

    return energy;
}

void get_LJ_DV(ATOM at[], DATA *dat, double fx[], double fy[], double fz[])
{
    uint32_t i=0 , j=0 ;
    double epsi_g=0.0 , sig_g=0.0;
    double dx=0.0 , dy=0.0 , dz=0.0 , d2=0.0 ;
    double de=0.0 ;

    for (i=0 ; i < dat->natom ; i++ )
    {
        fx[i] = 0.0 ;
        fy[i] = 0.0 ;
        fz[i] = 0.0 ;
        for (j=0 ; j < dat->natom ; j++ )
        {
            if (i==j) continue ;
            epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
            sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );
            dx = at[i].x - at[j].x ;
            dy = at[i].y - at[j].y ;
            dz = at[i].z - at[j].z ;
            d2  = dx*dx + dy*dy + dz*dz ;
            de = -24.0*epsi_g*( 2.0*(X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) ) / d2 ;
            fx[i] += de*dx;
            fy[i] += de*dy;
            fz[i] += de*dz;
        }
    }
}

double get_AZIZ_V(ATOM at[], DATA *dat, int32_t candidate)
{

    uint32_t i,j;
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
                        etmp=aziz_ne_ne(d);
                    else	//both are argon
                        etmp=aziz_ar_ar(d);

                }
                else //if different type it means it is 1 Ar and 1 Ne
                    etmp=aziz_ar_ne(d);

                energy += etmp;
            }
        }
    }
    else
    {
        i = (uint32_t) candidate;
        for (j=0; j<(dat->natom); j++)
        {
            if (j!=i)
            {
                d2 = X2(at[i].x-at[j].x) +  X2(at[i].y-at[j].y) + X2(at[i].z-at[j].z) ;
                d = sqrt(d2);

                if(!strcmp(at[i].sym,at[j].sym)) //if the same type (strcmp return 0 if identical)
                {
                    if (!strcmp(at[i].sym,"Ne")) //both are neon
                        etmp=aziz_ne_ne(d);
                    else	//both are argon
                        etmp=aziz_ar_ar(d);

                }
                else //if different type it means it is 1 Ar and 1 Ne
                    etmp=aziz_ar_ne(d);

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

//     return 0.0;

}

double aziz_ne_ne(double r)
{
    // HFD-B potential from Aziz, Chem. Phys. 130 (1989) p 187
    // reads r in angstroems gives potential energy in cm-1.

    // local
    double x;
    double v_star;
    double x2,x4,x6,x8,x10,f;

    // local constants
    const double a_s=8.9571795e+5;
    const double alpha_s=13.86434671;
    const double c_6=1.21317545;
    const double c_8=0.53222749;
    const double c_10=0.24570703;
    const double beta_s=-0.12993822;
    const double d=1.36;
    const double epsi=42.25*0.695039;
    const double r_m=3.091;

    // output
    double pot;


    x=r/r_m;

    f=exp(-X2(d/x-1.));

    if (x>=d)
        f=1.;

    x2=x*x;
    x4=x2*x2;
    x6=x2*x4;
    x8=x4*x4;
    x10=x4*x6;

    v_star=a_s*exp(-alpha_s*x+beta_s*X2(x))-f*(c_6/x6+c_8/x8+c_10/x10);

    pot=epsi*v_star;

    return pot;
}

double aziz_ar_ne(double r)
{
    // HFD-B potential from Barrow, Aziz, JCP 89, 6189 (1988)
    // reads r in angstroems gives potential energy in cm-1.
    
    // local
    double x;
    double v_star;
    double x2,x4,x6,x8,x10,x12,f;

    // local constants
    const double a_s=1.651205e+5;
    const double alpha_s=9.69290567;
    const double c_6=1.09781826;
    const double c_8=0.34284623;
    const double c_10=0.30103922;
    const double c_12=0.74483225;
    const double beta_s=-2.27380851;
    const double d=1.44;
    const double epsi=67.59*0.695039;
    const double r_m=3.4889;
    
    // received
    double pot;

    x=r/r_m;

    f=exp(-X2(d/x-1.));
    
    if (x>=d)
        f=1.;

    x2=x*x;
    x4=x2*x2;
    x6=x2*x4;
    x8=x4*x4;
    x10=x4*x6;
    x12=x6*x6;

    v_star=a_s*exp(-alpha_s*x+beta_s*X2(x))-f*(c_6/x6+c_8/x8+c_10/x10+c_12/x12);

    pot=epsi*v_star;

    return pot;
}

double aziz_ar_ar(double r)
{
    //HFD-B potential from Aziz, JCP 92, 1030 (1990)
    //reads r in angstroems gives potential energy in cm-1.
    
    // local
    double x;
    double v_star;
    double x2,x4,x6,x8,x10,f;

    // local constants
    const double a_s=1.14211845e+5;
    const double alpha_s=9.00053441;
    const double c_6=1.09971113;
    const double c_8=0.54511632;
    const double c_10=0.39278653;
    const double beta_s=-2.60270226;
    const double d=1.04;
    const double epsi=143.25*0.695039;
    const double r_m=3.761;
    
    // received
    double pot;

    x=r/r_m;

    f=exp(-X2(d/x-1.));
    
    if (x>=d)
        f=1.;

    x2=x*x;
    x4=x2*x2;
    x6=x2*x4;
    x8=x4*x4;
    x10=x4*x6;

    v_star=a_s*exp(-alpha_s*x+beta_s*X2(x))-f*(c_6/x6+c_8/x8+c_10/x10);

    pot=epsi*v_star;

    return pot;
}

