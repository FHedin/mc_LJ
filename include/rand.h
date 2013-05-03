#ifndef RAND_H_INCLUDED
#define RAND_H_INCLUDED

double get_next(DATA *dat);
double get_BoxMuller(DATA *dat, SPDAT *spdat);
void test_norm_distrib(DATA *dat, SPDAT *spdat, int n);

#endif // RAND_H_INCLUDED
