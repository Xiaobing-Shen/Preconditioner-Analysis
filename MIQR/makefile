#
OBJ= main.o auxf.o pcgnr.o matfun.o imqr.o tools.o
#
#LINKS = -L/project/darpa/lib/PC -llapack \
#        -L/project/darpa/lib/PC -lblas \
#	-L/project/darpa/lib/PC -lskit 
#
FC      =  f77
FCFLAGS =  -c -O2
CC      =  cc
CCFLAGS =  -c -DLINUX -O2
LD      =  f77 
LDFLAGS = 
LIB     = ${OBJblas} ${OBJeispack} ${OBJlapack}
#
# clear list of default suffixes, and declare default suffixes
.SUFFIXES:
.SUFFIXES: .f .c .o
# default rule to make .o files from .f files
.f.o  : ;       $(FC) $(FCFLAGS) $*.f -o $*.o
.c.o  : ;       $(CC) $(CCFLAGS) $*.c -o $*.o
#
mqr.ex: ${OBJ} ${LIB}
	$(LD)  $(LDFLAGS)  ${OBJ} ${LIB} $(LINKS) -o mqr.ex
#
clean :
	rm -f  ${OBJ} ${LIB} mqr.ex core 
