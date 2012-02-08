CC=gcc #-parallel
GSL_INC =# -I/opt/local/var/macports/software/gsl-devel/1.11.90_0/opt/local/include
GSL_LIB =# -L/opt/local/var/macports/software/gsl-devel/1.11.90_0/opt/local/lib
CFLAGS=-Wall -O2 $(OPT) $(GSL_INC)

#CC=icc -openmp
#CFLAGS= -Wall#-g -O2 -marchcore2 -w1 $(OPT)
#-marchcore2: optimized for Intel Dual Core

LIBS= -lm -lgsl -lgslcblas $(GSL_LIB)

EXEC = substructureOM#executables
OBJS = gadget2conv.o substructure.o#.o-files
INCL = Makefile gadget2conv.h#add headerfiles

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)  

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)

