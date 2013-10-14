#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

typedef struct
{
    char crdtitle_first[128];
    char crdtitle_last[128];
    
    char trajtitle[128];
    
    char etitle[128];
    
    uint32_t esave;
    uint32_t trsave;
} IODAT;

extern IODAT io;

extern FILE *crdfile;
extern FILE *traj;
extern FILE *efile;
extern FILE *stfile;

//pointer to the desired IO function
void (*write_traj)(ATOM at[], DATA *dat, uint64_t when);

void read_xyz(ATOM at[], DATA *dat, FILE *inpf);
void write_xyz(ATOM at[], DATA *dat, uint64_t when, FILE *outf);
void write_dcd(ATOM at[], DATA *dat, uint64_t when);

void write_rst(DATA *dat, SPDAT *spdat, ATOM at[]);

#endif // IO_H_INCLUDED
