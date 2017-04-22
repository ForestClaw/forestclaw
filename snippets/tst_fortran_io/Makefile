all:
	mpif90 -c hello.f90
	mpicc -g -c prog.c -lmpi
	mpicc -o tst_prog prog.o hello.o

clean:
	rm -f *.o *~ *.mod tst_prog
