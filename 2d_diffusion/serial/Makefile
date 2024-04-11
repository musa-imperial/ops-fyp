CC       = g++
OMPFLAGS = -Ofast -lm

BIN =  diffusion
LIBS = -lboost_chrono

all: $(BIN)

diffusion: diffusion.cpp Makefile
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ diffusion.cpp $(LIBS)

clean:
	$(RM) $(BIN)
	$(RM) *.txt

