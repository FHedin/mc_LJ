/*
 * Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 *
 */

#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

// structure containing files path and frequency of writting, for things related to the simulation (errors handled separately)
typedef struct
{
    // path for files where initial and final crds are stored in xyz format
    char crdtitle_first[FILENAME_MAX];
    char crdtitle_last[FILENAME_MAX];

    // path for file where the trajectory is stored
    char trajtitle[FILENAME_MAX];

    // path for file where the energy is stored
    char etitle[FILENAME_MAX];

    // frequency for saving energy and trajectory
    uint32_t esave;
    uint32_t trsave;
} IODAT;

// the previous structure is a global variable initialised in main.c
extern IODAT io;

// FILEs are opened in main.c and are global
extern FILE *crdfile;
extern FILE *traj;
extern FILE *efile;

//pointer to the desired IO function
void (*write_traj)(ATOM at[], DATA *dat, uint64_t when);

void read_xyz(ATOM at[], DATA *dat, FILE *inpf);
void write_xyz(ATOM at[], DATA *dat, uint64_t when, FILE *outf);
void write_dcd(ATOM at[], DATA *dat, uint64_t when);

void write_rst(DATA *dat, SPDAT *spdat, ATOM at[]);

#endif // IO_H_INCLUDED
