all:
	mpif90 -c hello.f90 hello_1.f90 -g
	mpicc -g -c prog.c -lmpi
	mpicc -o tst_prog prog.o hello.o hello_1.o

clean:
	rm -f *.o *~ *.mod tst_prog
