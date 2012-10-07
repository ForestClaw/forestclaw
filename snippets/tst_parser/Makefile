TARGET = parser_example

LDFLAGS = -L/Users/calhoun/projects/ForestClaw/code/p4est-build/local/lib -lsc
CFLAGS = -I.
CXXFLAGS = -I.
CC = c++

OBJS = example_parser.o dictionary.o parser.o

all: $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f *.o $(TARGET)
