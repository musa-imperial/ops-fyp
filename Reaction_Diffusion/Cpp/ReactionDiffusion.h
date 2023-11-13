#pragma once

class ReactionDiffusion {
private:
	// Private parameters
	double* u;
	double* v;
	double* u_calc; // array to store next values of u as they are calculated
	double* v_calc; // array to store next values of v as they are calculated

public:
	// Public parameters
	double dt, a, b, mu1, mu2, eps;
	int T, Nx, Ny;
	
	// Constructor
	ReactionDiffusion();

	// Destructor
	~ReactionDiffusion();
		
	// Setting parameters
	void SetParameters(const double& dt, const int& T, const int& Nx, const int& Ny, const double& a, const double& b, const double& mu1, const double& mu2, const double& eps);

	// Solving the problem
	void TimeIntegrate();
	
	// Setting initial conditions
	void SetInitialConditions();
	
	// Writing results to text file
	void WriteResults(const std::string& filename) const;
};
