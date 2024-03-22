CC       = g++
OMPFLAGS = -Ofast -lm

BIN =  diffusion
LIBS = 

all: $(BIN)

diffusion: diffusion.cpp Makefile
	$(CC) $(CCFLAGS) $(OMPFLAGS) -o $@ diffusion.cpp $(LIBS)

clean:
	$(RM) $(BIN)
	$(RM) *.txt

