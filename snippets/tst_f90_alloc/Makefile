all:
	gfortran -g -O2 -c array.f90 array_handling.f90
	gcc -g -O2 -c prog.c
	gcc -g -o tst_prog -O2 prog.o array.o array_handling.o -lgfortran

clean:
	rm -f *.o *~ *.mod tst_prog
