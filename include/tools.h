#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

void get_vector(DATA *dat,int mv_direction, double vec[3]);

void build_cluster(ATOM at[], DATA *dat, int from, int to, int mode);
int  no_conflict(ATOM at[],int i);

void steepd(ATOM at[],DATA *dat);

void adj_dmax(DATA *dat, int *step, int *acc);

CM getCM(ATOM at[],DATA *dat);
void recentre(ATOM at[], DATA *dat);

#endif // TOOLS_H_INCLUDED
