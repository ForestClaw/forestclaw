TARGET = parser_example

LDFLAGS = -L$(P4EST_DIR)/lib -lsc

CXXFLAGS = -I.
CXX = c++

CFLAGS = $(CXXFLAGS)
CC = gcc

OBJS = example_parser.o parser.o

all: $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f *.o $(TARGET)
