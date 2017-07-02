read:
	mpicc -c read_ncfile.c
	mpicc -o read_ncfile -O2 read_ncfile.o -lnetcdf -lhdf5_hl -lhdf5 -lcurl

write:
	mpicc -c write_ncfile.c
	mpicc -o write_ncfile -O2 write_ncfile.o -lnetcdf -lhdf5_hl -lhdf5 -lcurl

.PHONY: all read_ncfile write_ncfile

clean:
	rm -f *.o
	rm write_ncfile
	rm read_ncfile