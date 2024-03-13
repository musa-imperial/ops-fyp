#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
// #include <boost/chrono.hpp>
// #include <boost/program_options.hpp>
#include <time.h>
#include <omp.h>

int main(int argc, const char** argv)
{

  //Size along y
  int jmax = 4094;
  //Size along x
  int imax = 4094;
  //Size along x
  int iter_max = 100;

  int i, j;

  double pi  = 2.0 * asin(1.0);
  const double tol = 1.0e-6;
  double error     = 1.0;

  double *A;
  double *Anew;
  double *y0;

  A    = (double *)malloc((imax+2) * (jmax+2) * sizeof(double));
  Anew = (double *)malloc((imax+2) * (jmax+2) * sizeof(double));
  y0   = (double *)malloc((imax+2) * sizeof(double));

  memset(A, 0, (imax+2) * (jmax+2) * sizeof(double));

  // set boundary conditions
  // boost::chrono::high_resolution_clock::time_point t1 = boost::chrono::high_resolution_clock::now(); // start time
  clock_t begin = clock();
  for (int i = 0; i < imax+2; i++)
    A[(0)*(imax+2)+i]   = 0.0;

  for (int i = 0; i < imax+2; i++)
    A[(jmax+1)*(imax+2)+i] = 0.0;

  for (int j = 0; j < jmax+2; j++)
  {
    A[(j)*(imax+2)+0] = sin(pi * j / (jmax+1));
  }

  for (int j = 0; j < imax+2; j++)
  {
    A[(j)*(imax+2)+imax+1] = sin(pi * j / (jmax+1))*exp(-pi);
  }

  printf("Jacobi relaxation Calculation: %d x %d mesh\n", imax+2, jmax+2);

  int iter = 0;

  for (int i = 1; i < imax+2; i++)
    Anew[(0)*(imax+2)+i]   = 0.0;

  for (int i = 1; i < imax+2; i++)
    Anew[(jmax+1)*(imax+2)+i] = 0.0;

  for (int j = 1; j < jmax+2; j++)
    Anew[(j)*(imax+2)+0]   = sin(pi * j / (jmax+1));

  for (int j = 1; j < jmax+2; j++)
    Anew[(j)*(imax+2)+jmax+1] = sin(pi * j / (jmax+1))*expf(-pi);

  int nThreads;
  #pragma omp parallel default(shared)
	{
		nThreads = omp_get_num_threads();
   
  }
   printf("%d \n", nThreads);
#pragma omp parallel \
    default(none) shared(iter_max,jmax,imax, i, error, A, Anew, iter) \

{
  while ( error > tol && iter < iter_max )
  {
    error = 0.0;

    #pragma omp for private(j) schedule(static)
    for( j = 1; j < jmax+1; j++ )
    {
      for(i = 1; i < imax+1; i++)
      {
        Anew[(j)*(imax+2)+i] = 0.25f * ( A[(j)*(imax+2)+i+1] + A[(j)*(imax+2)+i-1]
            + A[(j-1)*(imax+2)+i] + A[(j+1)*(imax+2)+i]);
        error = fmax( error, fabs(Anew[(j)*(imax+2)+i]-A[(j)*(imax+2)+i]));
      }
    }
    #pragma omp for private(j) schedule(static)
    for(j = 1; j < jmax+1; j++ )
    {
      for(i = 1; i < imax+1; i++)
      {
        A[(j)*(imax+2)+i] = Anew[(j)*(imax+2)+i];    
      }
    }
    #pragma omp master 
    {
    if(iter % 10 == 0) printf("%5d, %0.6f\n", iter, error);   
    iter++;     
    }
    #pragma omp barrier
    
  }
}
  clock_t end = clock();

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("%5d, %0.6f\n", iter, error);

  // boost::chrono::high_resolution_clock::time_point t2 = boost::chrono::high_resolution_clock::now(); // end time
	
	std::cout << "-- Run-time: " << 
			time_spent << " ms --\n";

  double err_diff = fabs((100.0*(error/2.421354960840227e-03))-100.0);
  printf("Total error is within %3.15E %% of the expected error\n",err_diff);
  if(err_diff < 0.001)
    printf("This run is considered PASSED\n");
  else
    printf("This test is considered FAILED\n");


  free(A);
  free(Anew);
  return 0;
}

