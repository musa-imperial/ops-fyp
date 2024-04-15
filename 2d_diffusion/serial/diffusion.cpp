#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <boost/chrono.hpp>

#define IDX(I,J) ((J)*Nx + (I))

void printMatrix(const double* array, int rows, int cols, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Couldn't open the file " << filename << std::endl;
        return;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            outFile << array[i * cols + j] << " ";
        }
        outFile << std::endl;
    }
    outFile.close();
}

int main(int argc, const char** argv)
{
  double t;
  double dt   = 0.001;
  double T    = 1.0;
  
  int    Nx   = 40;
  int    Ny   = 40;
  int    Npts = Nx*Ny;

  double dx = 1;
  double dy = 1;
  double Lx   = dx*(Nx-1);
  double Ly   = dy*(Ny-1);

  // double Lx = 1;
  // double Ly = 1;
  // double dx = Lx / (Nx-1);
  // double dy = Ly / (Ny-1);

  double nu   = 0.1;

  //calculate dx, dy
  dx = Lx / (Nx-1);
  dy = Ly / (Ny-1);

  double hnudt = nu*dt/dx/dx;


  int i, j;

  double pi  = 2.0 * asin(1.0);

  double *A = nullptr;
  double *Anew = nullptr;

  double u, error, max_error;
  error = 0;
  max_error = 0;
  

  A   = new double[Npts];
  Anew   = new double[Npts];
  

  // set boundary conditions
  boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now(); // start time
  clock_t begin = clock();
  for (i = 0; i < Nx; i++)
    A[IDX(i,0)]   = 0;

  for (i = 0; i < Nx; i++)
    A[IDX(i,Ny-1)] = 0;

  for (j = 0; j < Ny; j++)
  {
    A[IDX(0,j)] = 0;//sin(pi * j / (Ny-1));
  }

  for (j = 0; j < Ny; j++)
  {
    A[IDX(Nx-1,j)] = 0;//sin(pi * j / (Ny-1))*exp(-pi);
  }

  for (i = 0; i < Nx; i++)
    Anew[IDX(i,0)]   = 0.0;

  for (i = 0; i < Nx; i++)
    Anew[IDX(i,Ny-1)] = 0.0;

  for (j = 0; j < Ny; j++)
  {
    Anew[IDX(0,j)] = sin(pi * j / (Ny-1));
  }

  for (j = 0; j < Ny; j++)
  {
    Anew[IDX(Nx-1,j)] = sin(pi * j / (Ny-1))*exp(-pi);
  }

  //initial condition

  for( j = 1; j < Ny-1; j++ )
    {
      for(i = 1; i < Nx-1; i++)
      {
        A[IDX(i,j)] = 5*sin(pi*(dx*i)/Lx)*sin(pi*(dy*j)/Ly);
        
        Anew[IDX(i,j)] = 5*sin(pi*(dx*i)/Lx)*sin(pi*(dy*j)/Ly);
        
      }
    }

  printMatrix(A, Ny, Nx, "initial_condition.txt");

  

  for (j = 0; j < Nx; j++)
  {
    Anew[IDX(Nx-1,j)] = sin(pi * j / (Ny-1))*exp(-pi);
  }

  for (t = 0; t < T; t+=dt)
  {
    for( j = 1; j < Ny-1; j++ )
    {
      for(i = 1; i < Nx-1; i++)
      {
        Anew[IDX(i,j)] = A[IDX(i,j)]+hnudt*(A[IDX(i+1,j)]+A[IDX(i-1,j)]+A[IDX(i,j+1)]+A[IDX(i,j-1)]-4*A[IDX(i,j)]); 
        
      }
    }
    
    for(j = 1; j < Ny-1; j++ )
    {
      for(i = 1; i < Nx-1; i++)
      {
        A[IDX(i,j)] = Anew[IDX(i,j)];    
      }
    }    
    
  }

  for( j = 1; j < Ny-1; j++ )
    {
      for(i = 1; i < Nx-1; i++)
      {
        u = 5*exp(-nu*pi*pi*(1/Lx/Lx+1/Ly/Ly)*T)*sin(pi/Lx*(dx*(i)))*sin(pi/Ly*(dy*(j)));

        error = error + sqrt(abs(A[IDX(i,j)]*A[IDX(i,j)]-u*u));
        max_error = fmax(max_error, fabs((A[IDX(i,j)]-u)/u));

        //std::cout << dx*(i) << ", "<< dy*(j) << ", " << (A[IDX(i,j)]) << ", " << u << ", " << fabs((A[IDX(i,j)]-u))/u << std::endl;
      }
    }

  boost::chrono::high_resolution_clock::time_point t2 = boost::chrono::high_resolution_clock::now();

  boost::chrono::milliseconds time_spent = boost::chrono::duration_cast<boost::chrono::milliseconds>(t2-t1);

  printMatrix(A, Ny, Nx, "output.txt");
	// std::cout << "-- Run-time: " << 
	// 		time_spent.count() << " ms --\n";

  std::cout << "\n" <<time_spent.count() << " " << max_error*100;

  delete[] A;
  delete[] Anew;
  return 0;
}

