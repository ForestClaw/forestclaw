TARGET = phase_iterative
FC = gfortran
F90 = gfortran
FFLAGS = -O3
LFLAGS = -L${HOME}/software/FISHPACK/fishpack4.1/lib -lfishpack
LINK = gfortran 


OBJECTS = \
	./phase_iterative.o \
	./cellave.o \
	./out2.o

.SUFFIXES: .o .f .f90


all: $(OBJECTS)
	$(LINK) $(OBJECTS) -o $(TARGET) $(LFLAGS)

phase_iterative.o: phase_iterative.f90
	$(F90) -c phase_iterative.f90

phase_cc.o: phase_cc.f
	$(F90) -c phase_cc.f

clean:
	rm *.o $(TARGET)

realclean :
	rm *.o $(TARGET) fort.*

### DO NOT remove this line - make depends on it ###
