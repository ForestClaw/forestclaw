all:
	gcc -g -c prog.c
	gfortran -g -c assign.f
	gcc -o tst_prog prog.o assign.o

cpp_tst:
	g++ -c -g prog.cpp
	gfortran -c -g assign.f
	g++ -o tst_prog prog.o assign.o

clean:
	rm -f *.o *~ *.mod tst_prog
