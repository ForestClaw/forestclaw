#
#		Makefile for the clawpack code:
#
#	To compile the claw1 library files type:  make lib
#       (from the AMRCLAWLIB sub-directory)
#
#	To make an executable, type:              make xclaw
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawhdf
#       (from the application sub-directory)
#
#	To compile a single file.f type:	  make file.o
#
#	To combine all source files type:	  make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
#	To combine all library source:    make backup
#       (date is appended)
#
FFLAGS = -O
LFLAGS =
F77 = gfortran -c
LINK	= gfortran

CLAW = ./

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

LIBOBJECTS = \
  $(CLAW)/claw1ez.o \
  $(CLAW)/bc1.o \
  $(CLAW)/setprob.o \
  $(CLAW)/setaux.o \
  $(CLAW)/b4step1.o \
  $(CLAW)/claw1.o \
  $(CLAW)/step1.o \
  $(CLAW)/step1fw.o \
  $(CLAW)/copyq1.o \
  $(CLAW)/inlinelimiter.o \
  $(CLAW)/limiter.o \
  $(CLAW)/philim.o \
  $(CLAW)/src1.o \
  $(CLAW)/rp1.o \
  $(CLAW)/out1.o

lib: $(LIBOBJECTS)

clean:
	-rm -f $(LIBOBJECTS)

### DO NOT remove this line - make depends on it ###
