#include "header.h"

ReactionDiffusion::ReactionDiffusion()
		: dt(0.001), a(0.75), b(0.06), mu1(5.0), mu2(0.0), eps(50.0), T(100), Nx(101), Ny(101)
	{
		std::cout << "-- ReactionDiffusion initialised with default parameters -- \n";
	}

ReactionDiffusion::~ReactionDiffusion() {
	// Deallocating all the heap-allocated arrays
	delete[] this->u;
	delete[] this->v;
	delete[] this->u_calc;
	delete[] this->v_calc;
}

void ReactionDiffusion::SetParameters(const double& dt, const int& T, const int& Nx, const int& Ny, const double& a, const double& b, const double& mu1, const double& mu2, const double& eps) {
	// Updating parameters with new values
	this->dt = dt;
	this->T = T;
	this->Nx = Nx;
	this->Ny = Ny;
	this->a = a;
	this->b = b;
	this->mu1 = mu1;
	this->mu2 = mu2;
	this->eps = eps;

	// Dislaying new parameter set
	std::cout << "\n-- New ReactionDiffusion parameters successfully set -- \n";
}

void ReactionDiffusion::SetInitialConditions() {
	// Allocate arrays on the heap
	this->u = new double[Nx*Ny];
	this->v = new double[Nx*Ny];
	this->u_calc = new double[Nx*Ny];
	this->v_calc = new double[Nx*Ny];
	
	std::cout << "-- Problem Parameters -- \n" <<
		"\n    dt: " << dt <<
		"\n     T: " << T <<
		"\n    Nx: " << Nx <<
		"\n    Ny: " << Ny <<
		"\n     a: " << a <<
		"\n     b: " << b <<
		"\n   mu1: " << mu1 <<
		"\n   mu2: " << mu2 <<
		"\n   eps: " << eps << "\n";

	std::cout << "\n-- Setting initial conditions -- \n";
	
	// Assuming dx = dy = 1
	//
	// |(x0, y0)   (x1, y0)  (x2, y0)  ...  (xNx, y0)|
	// |(x0, y1)   (x1, y1)  
	// |(x0, y2)
	// |...
	// |(x0, yNy)
	
	// Implementing u(x, y) = {1 if y > Ly/2; 0 otherwise} @ t = 0
	int Ly = Ny - 1;
	int yLim = Ly/2 + Ny%2;
	for (int i = 0; i < Nx; i++) {
		for (int j = yLim; j < Ny; j++) {
			u[i*Ny + j] = 1.0;
		}
	}
	
	// Implementing v(x, y) = {a/2 if x < Lx/2; 0 otherwise} @ t = 0
	int Lx = Nx - 1;
	int xLim = Lx/2 - Nx%2;
	for (int i = 0; i < xLim; i++) {
		for (int j = 0; j < Ny; j++) {
			v[i*Ny + j] = a * 0.5;
		}
	}
}

void ReactionDiffusion::WriteResults(const std::string& filename) const {
	std::cout << "-- Writing results to " << filename << " --\n";
	
	std::ofstream output_file(filename);
	
	if (!output_file.good()) {
		std::cout << "Error with opening file with name '" << filename << "'. Terminating program \n";
		exit(1);
	}
	
	while (output_file.is_open()) {
		for (int j = 0; j < Ny; ++j) {
			for (int i = 0; i < Nx; ++i) {
				output_file <<
					i << " " <<
					j << " " <<
					u[i*Ny + j] << " " <<
					v[i*Ny + j] << " \n";
			}
			output_file << "\n";
		}
		output_file.close();
	}
	
	std::cout << "-- Finished writing results. Integration end -- \n";
}

void ReactionDiffusion::TimeIntegrate() {
	std::cout << "\n-- Solving the problem up to time T = " << T << " with a time-step (dt) of " << dt << " and spacial step (h) of " << 1 << " --\n";
	
	// Tracking progress to print to terminal - disabled for submission
	// int prog = 0;
	
	// Precomputing end index of the grid array
	int end = Nx*Ny - 1;
	
	// Variables for the parallel code
	int i, j;
	int nThreads;
	double aInv = 1/a; // precomputing 1/a to avoid division
	int col;
	
	// Extracting number of threads one time before running
	#pragma omp parallel default(shared)
	{
		nThreads = omp_get_num_threads();
	}
	std::cout << "    Number of threads: " << nThreads << "\n";
	
	for (double t = 0.0; t < T; t += dt) {
		// x = 0, y = 0 node calculation for u and v
		u_calc[0] = u[0] + dt * (
			mu1 * (u[Ny] + u[1] - 2.0*u[0]) +
			eps*u[0] * (1.0 - u[0]) * (u[0] - (v[0] + b) * aInv));
		v_calc[0] = v[0] + dt * (
			mu2 * (v[Ny] + v[1] - 2.0*v[0]) +
			u[0]*u[0]*u[0] - v[0]);
		
		// x = 0 inner column calculation for u and v
		for (j = 1; j < Ny-1; ++j) {
			u_calc[j] = u[j] + dt * (
				mu1 * (u[j-1] + u[j+1] + u[Ny + j] - 3.0*u[j]) +
				eps*u[j] * (1.0 - u[j]) * (u[j] - (v[j] + b) * aInv));
		}
		for (j = 1; j < Ny-1; ++j) {
			v_calc[j] = v[j] + dt * (
				mu2 * (v[j-1] + v[j+1] + v[Ny + j] - 3.0*v[j]) +
				u[j]*u[j]*u[j] - v[j]);
		}
		
		// x = 0, y = Ny-1 node calculation for u and v
		u_calc[Ny-1] = u[Ny-1] + dt * (
			mu1 * (u[Ny-2] + u[2*Ny-1] - 2.0*u[Ny-1]) +
			eps*u[Ny-1] * (1.0 - u[Ny-1]) * (u[Ny-1] - (v[Ny-1] + b) * aInv));
		v_calc[Ny-1] = v[Ny-1] + dt * (
			mu2 * (v[Ny-2] + v[2*Ny-1] - 2.0*v[Ny-1]) +
			u[Ny-1]*u[Ny-1]*u[Ny-1] - v[Ny-1]);
		
		#pragma omp parallel \
			default(none) shared(i,aInv,end)
		{
			#pragma omp for private(j,col) schedule(static)
			for (i = 1; i < Nx-1; ++i) {
				col = i*Ny; // precomputing current column for use
				
				// y = 0 boundary calculation for current column
				u_calc[col] = u[col] + dt * (
					mu1 * (u[(i-1)*Ny] + u[col + 1] + u[(i+1)*Ny] - 3.0*u[col]) +
					eps*u[col] * (1.0 - u[col]) * (u[col] - (v[col] + b) * aInv));
				v_calc[col] = v[col] + dt * (
					mu2 * (v[(i-1)*Ny] + v[col + 1] + v[(i+1)*Ny] - 3.0*v[col]) +
					u[col]*u[col]*u[col] - v[col]);
				
				// Calculation of inside of the current column for u and v
				for (j = 1; j < Ny-1; ++j) {
					u_calc[col + j] = u[col + j] + dt * (
						mu1 * (u[(i-1)*Ny + j] + u[(i+1)*Ny + j] + u[col + j-1] + u[col + j+1] - 4.0*u[col + j]) +
						eps*u[col + j] * (1.0 - u[col + j]) * (u[col + j] - (v[col + j] + b) * aInv));
				}
				for (j = 1; j < Ny-1; ++j) {
					v_calc[col + j] = v[col + j] + dt * (
						mu2 * (v[(i-1)*Ny + j] + v[(i+1)*Ny + j] + v[col + j-1] + v[col + j+1] - 4.0*v[col + j]) +
						u[col + j]*u[col + j]*u[col + j] - v[col + j]);
				}
				
				// y = Ny-1 boundary calculation for current column
				u_calc[(i+1)*Ny - 1] = u[(i+1)*Ny - 1] + dt * (
					mu1 * (u[col - 1] + u[(i+1)*Ny - 2] + u[(i+2)*Ny - 1] - 3.0*u[(i+1)*Ny - 1]) +
					eps*u[(i+1)*Ny - 1] * (1.0 - u[(i+1)*Ny - 1]) * (u[(i+1)*Ny - 1] - (v[(i+1)*Ny - 1] + b) * aInv));
				v_calc[(i+1)*Ny - 1] = v[(i+1)*Ny - 1] + dt * (
					mu2 * (v[col - 1] + v[(i+1)*Ny - 2] + v[(i+2)*Ny - 1] - 3.0*v[(i+1)*Ny - 1]) +
					u[(i+1)*Ny - 1]*u[(i+1)*Ny - 1]*u[(i+1)*Ny - 1] - v[(i+1)*Ny - 1]);
			}
		}
		// END PARALLEL REGION
		
		// x = Nx-1, y = 0 node calculation for u and v
		u_calc[end-Ny + 1] = u[end-Ny + 1] + dt * (
			mu1 * (u[end-2*Ny + 1] + u[end-Ny + 2] - 2.0*u[end-Ny + 1]) +
			eps*u[end-Ny + 1] * (1.0 - u[end-Ny + 1]) * (u[end-Ny + 1] - (v[end-Ny + 1] + b) * aInv));
		v_calc[end-Ny + 1] = v[end-Ny + 1] + dt * (
			mu2 * (v[end-Ny + 1-Ny] + v[end-Ny + 1 + 1] - 2.0*v[end-Ny + 1]) +
			u[end-Ny + 1]*u[end-Ny + 1]*u[end-Ny + 1] - v[end-Ny + 1]);
		
		// x = Nx-1 inner column calculation for u and v
		for (j = 1; j < Ny-1; ++j) {
			u_calc[end-Ny + 1 + j] = u[end-Ny + 1 + j] + dt * (
				mu1 * (u[end-Ny + j] + u[end-Ny + j+2] + u[end-2*Ny + 1 + j] - 3.0*u[end-Ny + 1 + j]) +
				eps*u[end-Ny + 1 + j] * (1.0 - u[end-Ny + 1 + j]) * (u[end-Ny + 1 + j] - (v[end-Ny + 1 + j] + b) * aInv));
		}
		for (j = 1; j < Ny-1; ++j) {
			v_calc[end-Ny + 1 + j] = v[end-Ny + 1 + j] + dt * (
				mu2 * (v[end-Ny + j] + v[end-Ny + j+2] + v[end-2*Ny + 1 + j] - 3.0*v[end-Ny + 1 + j]) +
				u[end-Ny + 1 + j]*u[end-Ny + 1 + j]*u[end-Ny + 1 + j] - v[end-Ny + 1 + j]);
		}
		
		// x = Nx-1, y = Ny-1 node calculation for u and v
		u_calc[end] = u[end] + dt * (
			mu1 * (u[end-Ny] + u[end-1] - 2.0*u[end]) +
			eps*u[end] * (1.0 - u[end]) * (u[end] - (v[end] + b) * aInv));
		v_calc[end] = v[end] + dt * (
			mu2 * (v[end-Ny] + v[end-1] - 2.0*v[end]) +
			u[end]*u[end]*u[end] - v[end]);
		
		// Swapping arround array pointers 
		std::swap<double*>(u, u_calc);
		std::swap<double*>(v, v_calc);
		
		// Displaying progress - disabled for submission
		/*
		if ((int)(t/T * 100) != prog) {
			prog = (int)(t/T * 100);
			std::cout << prog << "%.. \n";
		}
		*/
	}
	
	std::cout << "-- Finished solving --\n";
}

