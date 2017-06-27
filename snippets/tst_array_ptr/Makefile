all:
	gfortran -O2 -c array.f90 -std=f2008
	gcc -O2 -c prog.c
	gcc -o tst_prog -O2 prog.o array.o -lgfortran

clean:
	rm -f *.o *~ *.mod tst_prog
