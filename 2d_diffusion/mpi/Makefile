CC       = mpicxx
OMPFLAGS = -O3 -lm

BIN =  diffusion
LIBS = -lboost_chrono

all: $(BIN)

diffusion: diffusion.cpp Makefile
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ diffusion.cpp $(LIBS)

test1: ./diffusion
		mpiexec -np 1 ./diffusion

test2: ./diffusion
		mpiexec -np 4 ./diffusion



clean:
	$(RM) $(BIN)
	$(RM) *.txt

