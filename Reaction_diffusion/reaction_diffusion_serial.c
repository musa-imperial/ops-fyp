//Original code developped by Joshua Binns
// auto-generated by ops.py
//
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
//#define OPS_2D
//#include <ops_seq_v2.h>

//#include "cw_kernels.h"

void printValues(double* u, int rows, int cols, char *filename) {

  FILE *file = fopen(filename, "w");

    // Check if the file is opened successfully
    if (file == NULL) {
        printf("Error: Unable to open the file.\n");
        //return 1; // Exit with an error code
    }
  
  for (int j = 0; j < rows; j++) {
    for (int i = 0; i < cols; i++) {
      fprintf(file, "%3f ", u[j * cols + i]);
    }
    fprintf(file, "\n");
  }

  fclose(file);
}


int main(int argc, const char** argv)
{
  
  
	// Tracking progress to print to terminal
	int prog = 0;
	
	// Precomputing end index of the grid array

  double dt = 0.01;
	double T = 200;
  printf("\n-- Solving the problem up to time T = %.2f with a time-step (dt) of %.2f and spacial step (h) of %d --\n", T, dt, 1);
	int    Nx = 101;
	int    Ny = 101;
	double a = 0.75;
	double b = 0.06;
	double mu1 = 5.0;
	double mu2 = 0.0;
	double eps = 50.0;

	int end = Nx*Ny - 1;
	
	// Variables for the parallel code
	int i, j;
	int nThreads;
	double aInv = 1/a; // precomputing 1/a to avoid division
	int col;
	
	double *u;
  double *u_calc;
  double *v;
  double *v_calc;

  u      = (double*)malloc((Nx*Ny) * sizeof(double));
  u_calc = (double*)malloc((Nx*Ny) * sizeof(double));
  v      = (double*)malloc((Nx*Ny) * sizeof(double));
  v_calc = (double*)malloc((Nx*Ny) * sizeof(double));
  memset(u, 0, (Nx*Ny)  * sizeof(double));
  memset(v, 0, (Nx*Ny)  * sizeof(double));


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

    for (int j = 0; j < Ny+1; j++) {
      for (int i = 0; i < Nx+1; i++) {
        u[j*Nx+i] = u_calc[j*Nx+i];
        v[j*Nx+i] = v_calc[j*Nx+i];
      }

    }
		
		// Displaying progress
		
		if ((int)(t/T * 100) != prog) {
			prog = (int)(t/T * 100);
      printf("%d%%...\n", prog);
      //printf("%.2f\n", t);
		}
		
	}
	
	printf("-- Finished solving --\n");

  printValues(u, Nx, Ny, "u_values");
  printValues(v, Nx, Ny, "v_values");

  return 0;
}

