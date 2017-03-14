target=prog

objs=prog.o global.o pkg.o

all:
	gcc -g -c prog.c global.c pkg.c -std=c99 -pedantic
	gcc -o ${target} ${objs}

clean:
	rm -f *.o *~ ${target} ${objs}
