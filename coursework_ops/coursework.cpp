#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, const char** argv)
{
  

  FILE *file = fopen("result.txt", "w");

  if (file == NULL) {
        printf("Unable to open the file.\n");
        return 1;
    }
//variables    
  double a = 1;
  double alpha = 0.003;
  double dx = 0.01;
  double dt = 0.001;
  double L = 1;

  int n = L/dx+1;
  double time = 0.06;
  double pi = 3.14159265359;

  double C1 = alpha*dt/dx/dx-a*dt/2/dx;
  double C2 = 1-2*alpha*dt/dx/dx;
  double C3 = a*dt/2/dx+alpha*dt/dx/dx;

  double *u;
  double *unew;

  u    = (double*)malloc((n) * sizeof(double));
  unew = (double*)malloc((n) * sizeof(double));

  memset(u, 0, (n)  * sizeof(double));

  // set boundary conditions
  u[0] = 0;

  for (int i = 1; i < n-1; i++) {
    
    for (int m=0; m < 5; m++)
        u[i] = u[i] + 10*pow(-1, m)*4/(2*m+1)/(2*m+1)/pi/pi*sin((2*m+1)*pi*(i*dx));

  }

  // int rows = 1;
  // int cols = n;
  // for (int i = 0; i < rows; i++) {
  //       for (int j = 0; j < n; j++) {
  //           fprintf(file, "%lf ", u[i * cols + j]);
  //       }
  //       fprintf(file, "\n");
  //   }
    
  
  
  u[n-1] = u[n-2];

  //time stepping 
  double t = 0;

  while ( t < time )
  {
    
    for( int i = 1; i < n-1; i++ )
    {
      unew[i] = C1*u[i+1]+C2*u[i]+C3*u[i-1];
    }
    
    for (int i = 1; i < n-1; i++) {
            u[i] = unew[i];
        }
    u[n-1] = u[n-2];
  // int rows = 1;
  // int cols = n;
  // for (int i = 0; i < rows; i++) {
  //       for (int j = 0; j < n; j++) {
  //           fprintf(file, "%lf ", u[i * cols + j]);
  //       }
  //       fprintf(file, "\n");
  // }
    t = t + dt;
  }

  int rows = 1;
  int cols = n;
  for (int i = 0; i < rows; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%lf ", u[i * cols + j]);
        }
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);

    printf("Matrix data has been written to result.txt.\n");





  printf("Program complete");

  

  free(u);
  free(unew);
  return 0;
}

