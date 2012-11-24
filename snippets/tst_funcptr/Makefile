Target := tst_funcptr

# default: $(TARGETS)

F77 = gfortran
CXX = g++
CFLAGS = -Wall -O0 -g

include Makefile.inc

OBJS =  $(Target).o \
	math_value.o \
	$(set_file) \
	$(not_defined) \
	$(math_file)

all: $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o $(Target)

.f.o:
	$(F77) $(CFLAGS) $^ -c -o $@

.cpp.o:
	$(CXX) $(CFLAGS) $^ -c -o $@

.c.o:
	$(CXX) $(CFLAGS) $^ -c -o $@

clean:
	rm -f *.o $(Target)

.PHONY: clean default
