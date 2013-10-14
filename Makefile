#################################################################
########################   MakeVars   ###########################
#################################################################

CC=gcc 
FC=gfortran
 
CC_OPT=-I"./dSFMT" -I"./include" -Wall -Wextra -std=c99 -O2 -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=19937

CC_SFMT_OPT=-I"./dSFMT" -Wall -Wextra -std=c99 -O2 -msse2 -fno-strict-aliasing -DHAVE_SSE2 -DDSFMT_MEXP=19937

FC_OPT=-O2 -msse2

LD_OPT=-lm

MKDIR=mkdir -p ./obj/dSFMT
 
CIBLE=myMC
 
SRC=$(wildcard ./src/*.c)
dSRC=$(wildcard ./dSFMT/*.c)
fSRC=$(wildcard ./src/*.f)
 
OBJ=$(patsubst ./src/%.c,./obj/%.o,$(SRC))
dOBJ=$(patsubst ./dSFMT/%.c,./obj/dSFMT/%.o,$(dSRC)) 
fOBJ=$(patsubst ./src/%.f,./obj/%.o,$(fSRC))
 
#################################################################
########################   Makefile   ###########################
#################################################################
 
all:$(CIBLE)
	@echo "Compilation Success"

$(CIBLE):Makefile

./obj/%.o:./src/%.c 
	$(CC) $(CC_OPT) -c $< -o $@ 

./obj/%.o:./src/%.f 
	$(FC) $(FC_OPT) -c $< -o $@ 

./obj/dSFMT/%.o:./dSFMT/%.c
	@$(MKDIR)
	$(CC) $(CC_SFMT_OPT) -c $< -o $@
 
$(CIBLE):$(dOBJ) $(fOBJ) $(OBJ)
	$(CC) $(dOBJ) $(fOBJ) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(CIBLE) ./obj/*.o

clean_all:
	rm -f $(CIBLE) ./obj/*.o ./obj/dSFMT/*.o

