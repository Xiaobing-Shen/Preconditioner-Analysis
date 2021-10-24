#
OBJ= main.o auxf.o pcgnr.o matfun.o imqr.o tools.o
#
BIT     =  64
FC      =  xlf
FCFLAGS =  -c -w -O2 -q$(BIT)
CC      =  cc
CCFLAGS =  -c -w -D_IBM -O2 -q$(BIT)
LD      =  $(FC)
LDFLAGS =  -q$(BIT)
LIB     = 
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
