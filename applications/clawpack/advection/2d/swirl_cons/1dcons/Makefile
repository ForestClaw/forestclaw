#
#               Makefile for the lib code:
#
#       To make an executable, type:            make xclaw
#       (from the application sub-directory)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To compile the lib library files type:    make claw1
#       (from the claw/clawpack/1d/lib sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)

#
#
FFLAGS =  -g
LFLAGS =  -O
LINK	= gfortran

CLAW = .

OBJECTS = \
  driver.o \
  qinit.o \
  rp1fwave.o \
  setaux.o \
  setprob.o

LIBOBJECTS = \
  ./lib/claw1ez.o \
  ./lib/b4step1.o \
  ./lib/bc1.o \
  ./lib/out1.o \
  ./lib/claw1.o \
  ./lib/step1fw.o \
  ./lib/src1.o \
  ./lib/copyq1.o \
  ./lib/limiter.o \
  ./lib/philim.o

xclaw: $(OBJECTS) $(LIBOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) -o xclaw

clean:
	rm -f *.o

### DO NOT remove this line - make depends on it ###
