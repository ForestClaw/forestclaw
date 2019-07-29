all:
	gfortran -g -c array.f
	gcc -g -c prog.c
	gcc -o tst_prog  prog.o array.o

clean:
	rm -f *.o *~ tst_prog
