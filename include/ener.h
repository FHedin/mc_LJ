#ifndef ENER_H_INCLUDED
#define ENER_H_INCLUDED

/*
 * Lennard Jones : divide by this for transforming
 * R_min in sigma ( it is 2^(1/6) )
 * */
#define RTOSIG 1.122462048309373

/*
 * Boltzmann cst in charmm units
 */
#define KBCH 1.98719e-03

/*
 * for converting energy
 */
#define JTOCAL      0.239005736     // Joules to Calories
#define CM1TOKJM    1.1963e-02    // 1 cm-1 in kJ/mol

//pointers to the desired energy and force functions
double (*get_ENER)(ATOM at[], DATA *dat, int candidate);
void   (*get_DV)(ATOM at[], DATA *dat, double DV[]);

//ener and force for lennard-jones
double get_LJ_V(ATOM at[], DATA *dat, int candidate);
void   get_LJ_DV(ATOM at[], DATA *dat, double DV[]);

//ener for aziz potential
double get_AZIZ_V(ATOM at[], DATA *dat, int candidate);
//those 3 fortran subroutines returns energy in cm-1
extern void aziz_ne_ne_(double *r, double *pot);
extern void aziz_ar_ne_(double *r, double *pot);
extern void aziz_ar_ar_(double *r, double *pot);

// constraint for avoiding cluster evaporation
double getExtraPot(double d2, double sig, double eps);

#endif // ENER_H_INCLUDED
