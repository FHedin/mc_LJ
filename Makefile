#
# Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
# All rights reserved.
#
# The 3-clause BSD license is applied to this software.
# see LICENSE.txt
# 
#

#################################################################
########################   MakeVars   ###########################
#################################################################

# intel compiler
CC=icc 
FC=ifort

#CC=gcc
#FC=gfortran

# CC=gcc48
# FC=gfortran48

# For compiling MS-Windows executable on linux 
# CC=x86_64-w64-mingw32-gcc
# FC=x86_64-w64-mingw32-gfortran

WFLAGS=-Wall -Wextra 
#-Wdouble-promotion -Wformat -Wimplicit-int -Wuninitialized -Wfloat-equal \
#-Wpointer-arith -Wtype-limits -Wbad-function-cast -Wcast-qual -Wconversion \
#-Wsign-conversion

OPTIM=-O2

INC_OPT=-I"./dSFMT/" -I"./include/" -I"/home/hedin/bin/luajit_last/include/luajit-2.0"
CC_OPT=$(INC_OPT) $(WFLAGS) -std=c99 $(OPTIM) -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=19937 -DLUA_PLUGINS

CC_SFMT_OPT=-I"./dSFMT" $(WFLAGS) -std=c99 $(OPTIM) -msse2 -fno-strict-aliasing -DHAVE_SSE2 -DDSFMT_MEXP=19937

FC_OPT=$(WFLAGS) -std=f95 $(OPTIM) -msse2

LD_OPT= -L"/home/hedin/bin/luajit_last/lib" -lluajit-5.1 -lm
#LD_OPT= -L"." -llua51 -lm
#LD_OPT= -lm

MKDIR=mkdir -p ./obj/dSFMT
 
TARGET=mc_LJ
 
SRC=$(wildcard ./src/*.c)
dSRC=$(wildcard ./dSFMT/*.c)
fSRC=$(wildcard ./src/*.f90)
 
OBJ=$(patsubst ./src/%.c,./obj/%.o,$(SRC))
dOBJ=$(patsubst ./dSFMT/%.c,./obj/dSFMT/%.o,$(dSRC)) 
fOBJ=$(patsubst ./src/%.f90,./obj/%.o,$(fSRC))
 
#################################################################
########################   Makefile   ###########################
#################################################################
 
all:$(TARGET)
	@echo "Compilation Success"

$(TARGET):Makefile

./obj/%.o:./src/%.c 
	$(CC) $(CC_OPT) -c $< -o $@ 

./obj/%.o:./src/%.f90 
	$(FC) $(FC_OPT) -c $< -o $@ 

./obj/dSFMT/%.o:./dSFMT/%.c
	@$(MKDIR)
	$(CC) $(CC_SFMT_OPT) -c $< -o $@
 
$(TARGET):$(dOBJ) $(fOBJ) $(OBJ)
	$(CC) $(CC_OPT) $(dOBJ) $(fOBJ) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(TARGET) ./obj/*.o

clean_all:
	rm -f $(TARGET) ./obj/*.o ./obj/dSFMT/*.o

