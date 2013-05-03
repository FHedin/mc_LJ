#ifndef MCCLASSIC_H_INCLUDED
#define MCCLASSIC_H_INCLUDED

int make_MC_moves(ATOM at[], DATA *dat, double *ener);
int apply_Metrop(ATOM at[], ATOM at_new[], DATA *dat, int *candidate, double *ener);

#endif // MCCLASSIC_H_INCLUDED
