CC       = g++
OMPFLAGS = -O3 -fopenmp 

BIN =  diffusion
LIBS = -lboost_chrono -lm

all: $(BIN)

diffusion: diffusion.cpp Makefile
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ diffusion.cpp $(LIBS)

clean:
	$(RM) $(BIN)
	$(RM) *.txt

