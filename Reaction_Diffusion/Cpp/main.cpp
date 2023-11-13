#include "header.h"

int main(int argc, char* argv[]) {
	boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now(); // start time
	
	// Initialise system
	ReactionDiffusion System;
	
	// Change parameters as dictated by command line arguments
	input_handler(argc, argv, System);
	
	System.SetInitialConditions();
	
	System.TimeIntegrate();
	
	System.WriteResults("output.txt");
	
	boost::chrono::high_resolution_clock::time_point t2 = boost::chrono::high_resolution_clock::now(); // end time
	
	std::cout << "-- Run-time for TimeIntegrate function: " << 
			boost::chrono::duration_cast<boost::chrono::milliseconds>(t2-t1) << " ms --\n";
}
