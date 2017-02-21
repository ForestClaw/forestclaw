all:
	gfortran -g -c array.f90 -std=f95
	gcc -g -c prog.c
	gcc -o tst_prog prog.o array.o -lgfortran

clean:
	rm -f *.o *~ *.mod tst_prog
