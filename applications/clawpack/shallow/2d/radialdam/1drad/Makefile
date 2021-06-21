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

F77 = gfortran
FC = gfortran
FFLAGS =  -O
LFLAGS =  -O
LINK	= gfortran

CLAW = .

OBJECTS = \
  driver.o \
  qinit.o \
  rp1sw.o \
  src1.o \
  setprob.o

LIBOBJECTS = \
  $(CLAW)/lib/setaux.o \
  $(CLAW)/lib/claw1ez.o \
  $(CLAW)/lib/b4step1.o \
  $(CLAW)/lib/bc1.o \
  $(CLAW)/lib/out1.o \
  $(CLAW)/lib/claw1.o \
  $(CLAW)/lib/step1.o \
  $(CLAW)/lib/copyq1.o \
  $(CLAW)/lib/limiter.o \
  $(CLAW)/lib/philim.o

xclaw: $(OBJECTS) $(LIBOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) -o xclaw

### DO NOT remove this line - make depends on it ###
