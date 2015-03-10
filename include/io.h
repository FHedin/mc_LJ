/**
 * \file io.h
 * 
 * \brief Header file : I/O related functions, structures and a few global variables
 * 
 * \copyright Copyright (c) 2011-2015, Florent Hédin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#ifndef IO_H_INCLUDED
#define IO_H_INCLUDED

#ifndef FILENAME_MAX
#define FILENAME_MAX    4096
#endif

/// structure containing files path and frequency of writing, for things related to the simulation (errors handled separately)
typedef struct
{
    /// path for file where initial crds are stored in xyz format
    char crdtitle_first[FILENAME_MAX];
    
    /// path for file final crds are stored in xyz format
    char crdtitle_last[FILENAME_MAX];

    /// path for file where the trajectory is stored
    char trajtitle[FILENAME_MAX];

    /// path for file where the energy is stored
    char etitle[FILENAME_MAX];

    /// frequency for saving energy
    uint32_t esave;
    /// frequency for saving trajectory
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

//read or write coordinates or trajectory files
void read_xyz(ATOM at[], DATA *dat, FILE *inpf);
void write_xyz(ATOM at[], DATA *dat, uint64_t when, FILE *outf);
void write_dcd(ATOM at[], DATA *dat, uint64_t when);

// BUG : restart file 
void write_rst(ATOM at[], DATA *dat, SPDAT *spdat, uint32_t meth);

#endif // IO_H_INCLUDED
