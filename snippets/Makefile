
TARGETS = tst_funcptr

default: $(TARGETS)

F77 = gfortran
CXX = g++
CFLAGS = -Wall -O0 -g

tst_funcptr: tst_funcptr.o set_value.o
	$(CXX) $(CFLAGS) $^ -o $@

.f.o:
	$(F77) $(CFLAGS) $^ -c -o $@

.cpp.o:
	$(CXX) $(CFLAGS) $^ -c -o $@

clean:
	rm -f *.a *.o $(TARGETS)

.PHONY: clean default
