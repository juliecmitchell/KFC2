###############################################################################
# Platform-dependent parameters
CC	= cc
CPP	= c++
F77	= f77
MPIC	=
MPIF	= hf77
CFLAGS	= -g -O2 
FFLAGS	= -O2
LDFLAGS	=  
LIBS	= -lm
SYS	=
SPECIAL	= 
RANLIB	= /bin/true
SHELL	= /bin/sh
###############################################################################

TARGETS	=	core_rim_1 core_rim_2 mut_compensation_1 splitpdb

CR1 	=	core_rim_1.o
CR2 	=	core_rim_2.o
MC1 	=	mut_compensation_1.o
SPDB 	=	splitpdb.o

MYINC	= 	core_rim_1.h core_rim_2.h mut_comp.h

#: Default rule
default: install clean   

# Rules to link object files into executables

core_rim_1: core_rim_1.o core_rim_1.h
	$(CC) -o core_rim_1.c $(LDFLAGS) $(LIBS)

core_rim_2: $(CR2) core_rim_2.h
	$(CC) -o $@ $(CR2) $(LDFLAGS) $(LIBS)

mut_compensation_1: $(MC1) mut_comp.h
	$(CC) -o $@ $(MC1) $(LDFLAGS) $(LIBS)

splitpdb: $(SPDB)
	$(CC) -o $@ $(SPDB) $(LDFLAGS) $(LIBS)

# Rule to build object file.o from file.c
.c.o:
	$(CC) -c $(CFLAGS) $(MYINC) $<

# Rule to build object file.o from file.C
.C.o:
	$(CPP) -c $(FFLAGS) $(MYINC) $<

# Rule to build object file.o from file.f
.f.o:
	$(F77) -c $(FFLAGS) $<

install:
	for executable in $(TARGETS); do \
		cc -o $${executable}.c $${executable}; \
		mv $${executable} ../bin/KFC2_exec/$${executable}; \
		chmod 755 ../bin/KFC2_exec/$${executable}; \
	done

clean:
	/bin/rm -f *.o *~
	/bin/rm -f $(TARGETS)

realclean: clean
	for executable in $(TARGETS); do \
		rm -f ../bin/KFC2_exec/$${executable}*; \
	done
