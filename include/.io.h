#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

extern IODAT io;

extern FILE *traj;
extern FILE *efile;
extern FILE *stfile;

extern int dcd_header_empty;

//pointer to the desired IO function
void (*write_traj)(ATOM at[], DATA *dat, int when);

void write_xyz(ATOM at[], DATA *dat, int when);
void write_dcd(ATOM at[], DATA *dat, int when);

#endif // IO_H_INCLUDED
