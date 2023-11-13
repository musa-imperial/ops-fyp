#include "header.h"

namespace po = boost::program_options;

void input_handler(int& argc, char** argv, ReactionDiffusion& System) {
			
		// Declare the supported options
		po::options_description desc("Allowed options");
		desc.add_options()
			("help", "display help message")
			("dt", po::value<double>(), "set time step")
			("T", po::value<int>(), "set final time")
			("Nx", po::value<int>(), "set number of x grid points")
			("Ny", po::value<int>(), "set number of y grid points")
			("a", po::value<double>(), "set value of a")
			("b", po::value<double>(), "set value of b")
			("mu1", po::value<double>(), "set value of mu1")
			("mu2", po::value<double>(), "set value of mu2")
			("eps", po::value<double>(), "set value of epsilon")
		;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
		
		// Current values stored in the System
		double dt = System.dt;
		double T = System.T;
		double Nx = System.Nx;
		double Ny = System.Ny;
		double a = System.a;
		double b = System.b;
		double mu1 = System.mu1;
		double mu2 = System.mu2;
		double eps = System.eps;
		
		// Tracking if any variables are changed
		bool isChanged = false;

		if (vm.count("help")) {
			std::cout << desc << "\n";
			exit(1);
		}
		
		// Check to see if any relevant command line arguments are used and update values
		if (vm.count("dt")) {
			dt = vm["dt"].as<double>();
			isChanged = true;
		}
		if (vm.count("T")) {
			T = vm["T"].as<int>();
			isChanged = true;
		}
		if (vm.count("Nx")) {
			Nx = vm["Nx"].as<int>();
			isChanged = true;
		}
		if (vm.count("Ny")) {
			Ny = vm["Ny"].as<int>();
			isChanged = true;
		}
		if (vm.count("a")) {
			a = vm["a"].as<double>();
			isChanged = true;
		}
		if (vm.count("b")) {
			b = vm["b"].as<double>();
			isChanged = true;
		}
		if (vm.count("mu1")) {
			mu1 = vm["mu1"].as<double>();
			isChanged = true;
		}
		if (vm.count("mu2")) {
			mu2 = vm["mu2"].as<double>();
			isChanged = true;
		}
		if (vm.count("eps")) {
			eps = vm["eps"].as<double>();
			isChanged = true;
		}
		
		// Update instance of ReactionDiffusion if parameters are updated
		if (isChanged) System.SetParameters(dt, T, Nx, Ny, a, b, mu1, mu2, eps);
}
