simple_xy_wr:
	gcc-6 -c simple_xy_wr.c
	gcc-6 -o simple_xy_wr -O2 simple_xy_wr.o -lnetcdf

all:
	gcc-6 -c main.c
	gcc-6 -o main -O2 main.o -lnetcdf

clean:
	rm -f *.o
	rm simple_xy_wr
	rm main