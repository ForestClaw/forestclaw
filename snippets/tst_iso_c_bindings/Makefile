all:
	gcc -c prog.c
	gfortran -c assign.f
	gcc -o tst_prog prog.o assign.o

cpp_tst:
	g++ -c prog.cpp
	gfortran -c assign.f
	g++ -o tst_prog prog.o assign.o

clean:
	rm -f *.o *~ tst_prog
