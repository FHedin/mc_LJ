#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

void get_vector(DATA *dat,int32_t mv_direction, double vec[3]);

void build_cluster(ATOM at[], DATA *dat, uint32_t from, uint32_t to, int32_t mode);
int32_t  no_conflict(ATOM at[],uint32_t i);

void steepd(ATOM at[],DATA *dat);
void steepd_ini(ATOM at[],DATA *dat);

void adj_dmax(DATA *dat, uint64_t *step, uint64_t *acc);

CM getCM(ATOM at[],DATA *dat);
void recentre(ATOM at[], DATA *dat);

#endif // TOOLS_H_INCLUDED
