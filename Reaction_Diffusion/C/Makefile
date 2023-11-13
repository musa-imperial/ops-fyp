CXX = g++
CXXFLAGS = -fopenmp -Wall -O3
HDRS = header.h ReactionDiffusion.h
OBJS = main.o input-handler.o ReactionDiffusion.o
LDLIBS = -lboost_chrono -lboost_program_options

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

BarkleyModel: $(OBJS)
	$(CXX) -fopenmp -o $@ $^ $(LDLIBS)

all: BarkleyModel

test1: ./BarkleyModel
	OMP_NUM_THREADS=1 ./BarkleyModel --Nx 101 --Ny 101 --a 0.75 --b 0.06 --eps 50.0 --mu1 5.0 --mu2 0.0

test2: ./BarkleyModel
	OMP_NUM_THREADS=1 ./BarkleyModel --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0

test3: ./BarkleyModel
	OMP_NUM_THREADS=1 ./BarkleyModel --Nx 101 --Ny 101 --a 0.5 --b 0.1 --eps 50.0 --mu1 5.0 --mu2 0.0

test4: ./BarkleyModel
	OMP_NUM_THREADS=1 ./BarkleyModel --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01

RunTime-b: ./BarkleyModel
	for nTests in 1 2 3 4 5 6 7 8 ; do \
		OMP_NUM_THREADS=$(n) ./BarkleyModel --Nx 251 --Ny 251 --a 0.75 --b 0.06 --eps 13.0 --mu1 5.0 --mu2 0.0	; \
		sleep 2 ; \ # pauses execution to mitigate prior run affecting next run
	done

RunTime-d: ./BarkleyModel
	for nTests in 1 2 3 4 5 6 7 8 ; do \
		OMP_NUM_THREADS=$(n) ./BarkleyModel --Nx 151 --Ny 81 --a 0.75 --b 0.0001 --eps 12.5 --mu1 1.0 --mu2 0.01 ; \
		sleep 2 ; \ # pauses execution to mitigate prior run affecting next run
	done

.PHONY: clean
	target

clean:
	-rm -f *.o BarkleyModel
