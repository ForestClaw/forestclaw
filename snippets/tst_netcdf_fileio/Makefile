all:
	mpicc -c main.c
	mpicc -o main -O2 main.o -lnetcdf -lhdf5_hl -lhdf5 -lcurl

clean:
	rm -f *.o
	rm main