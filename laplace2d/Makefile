CC       = g++
OMPFLAGS = -Ofast -lm

BIN =  laplace2d
LIBS = -fopenmp #-lboost_chrono -lboost_program_options

all: $(BIN)

laplace2d: laplace2d_original.cpp Makefile
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ laplace2d_original.cpp $(LIBS)

clean:
	$(RM) $(BIN)

