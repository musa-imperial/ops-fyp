
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define OPS_2D
#include <ops_seq_v2.h>


int imax, jmax;
double pi  = 2.0 * asin(1.0);
#include "laplace_kernels.h"
// void ops_par_loop_set_zero(char const *, ops_block, int , int*,
//   ops_arg );

// void ops_par_loop_left_bndcon(char const *, ops_block, int , int*,
//   ops_arg,
//   ops_arg );

// void ops_par_loop_right_bndcon(char const *, ops_block, int , int*,
//   ops_arg,
//   ops_arg );

// void ops_par_loop_apply_stencil(char const *, ops_block, int , int*,
//   ops_arg,
//   ops_arg,
//   ops_arg );

// void ops_par_loop_copy(char const *, ops_block, int , int*,
//   ops_arg,
//   ops_arg );

int main(int argc, const char** argv)
{
  //Initialise the OPS library, passing runtime args, and settings diagnostics level to low (1)
  ops_init(argc, argv, 1);
  
  //Size along y
  jmax = 510;//4094;
  //Size along x
  imax = 510;//4094;
  //Size along x
  int iter_max = 100;

  // FILE *fptr = fopen("/home/musa/apps/OPS/apps/c/laplace2d_tutorial/tutorial/results.txt", "w");

  // if(fptr = NULL) {
  //   printf("Error!");
  //   exit(1);
  // }

  // const char *filename = "results.txt";
  // FILE *file = fopen(filename, "w"); // Open the file in write mode ("w")


  // if (file == NULL) {
  //       perror("Error opening file");
  //       return 1;
  //   }
  
  // fprintf(file, "Hello, World\n");

  
  //float pi  = 2.0 * asin(1.0);
  const double tol = 1.0e-10;
  double error     = 1.0;

  double *A=NULL;
  double *Anew=NULL;
  //_Float16 *y0;

  // A    = (_Float16 *)malloc((imax+2) * (jmax+2) * sizeof(_Float16));
  // Anew = (_Float16 *)malloc((imax+2) * (jmax+2) * sizeof(_Float16));
  // y0   = (_Float16 *)malloc((imax+2) * sizeof(_Float16));

  //memset(A, 0, (imax+2) * (jmax+2) * sizeof(_Float16));

  //The 2D block
  ops_block block = ops_decl_block(2, "my_grid");

  //The two datasets

  int size[] = {imax, jmax}; //size of grid or mesh points in each dimension, c = {0, 0}, fortran = {1, 1}

  int base[] = {0,0}; //think of it as what index is the starting index
  int d_m[] = {-1,-1}; //negative x and y direction for boundaries ?
  int d_p[] = {1,1}; //positive x and y directions for boundaries ?


  ops_dat d_A = ops_decl_dat(block, 1, size, base,
                              d_m, d_p, A, "double", "A");

  ops_dat d_Anew = ops_decl_dat(block, 1, size, base, 
                              d_m, d_p, Anew, "double", "Anew");

  ops_reduction h_err = ops_decl_reduction_handle(sizeof(double), "double", "error");
  
  //declar and define global constants

  ops_decl_const("imax", 1, "int", &imax);
  ops_decl_const("jmax", 1, "int", &jmax);
  ops_decl_const("pi", 1, "double", &pi);

  //Two stencils, a 1-points, and a 5-point

  int s2d_00[] = {0, 0};
  ops_stencil S2D_00 = ops_decl_stencil(2, 1, s2d_00, "0, 0");

  int s2d_5pt[] = {0,0, 1,0, -1,0, 0,1, 0,-1};
  ops_stencil S2D_5pt = ops_decl_stencil(2, 5, s2d_5pt, "5pt");

  ops_partition("");
  // set boundary conditions
  //for (int i = 0; i < imax+2; i++)
  //  A[(0)*(imax+2)+i]   = 0.0;

  //bottom bc
  int bottom_range[] = {-1, imax+1, -1, 0};

  int top_range[] = {-1, imax+1, jmax, jmax+1};

  int left_range[] = {-1, 0, -1, jmax+1};

  int right_range[] = {imax, imax+1, -1, jmax+1};

  int interior_range[] = {0, imax, 0, jmax};


  ops_par_loop(set_zero, "set_zero", block, 2, bottom_range,
                ops_arg_dat(d_A, 1, S2D_00, "double", OPS_WRITE));

  ops_par_loop(set_zero, "set_zero", block, 2, top_range,
                ops_arg_dat(d_A, 1, S2D_00, "double", OPS_WRITE));

  ops_par_loop(left_bndcon, "left_bndcon", block, 2, left_range,
              ops_arg_dat(d_A, 1, S2D_00, "double", OPS_WRITE),
              ops_arg_idx());

  ops_par_loop(right_bndcon, "right_bndcon", block, 2, right_range,
              ops_arg_dat(d_A, 1, S2D_00, "double", OPS_WRITE),
              ops_arg_idx());


  ops_par_loop(set_zero, "set_zero", block, 2, bottom_range,
                ops_arg_dat(d_Anew, 1, S2D_00, "double", OPS_WRITE));

  ops_par_loop(set_zero, "set_zero", block, 2, top_range,
                  ops_arg_dat(d_Anew, 1, S2D_00, "double", OPS_WRITE));


  ops_par_loop(left_bndcon, "left_bndcon", block, 2, left_range,
              ops_arg_dat(d_Anew, 1, S2D_00, "double", OPS_WRITE),
              ops_arg_idx());

  ops_par_loop(right_bndcon, "right_bndcon", block, 2, right_range,
              ops_arg_dat(d_Anew, 1, S2D_00, "double", OPS_WRITE),
              ops_arg_idx());

  // ops_par_loop(copy, "copy", block, 2, interior_range,
  //             ops_arg_dat(d_A,    1, S2D_00, "double", OPS_WRITE),
  //             ops_arg_dat(d_Anew, 1, S2D_00, "double", OPS_READ));
  

  // ops_par_loop(apply_stencil, "apply_stencil", block, 2, interior_range, 
  //         ops_arg_dat(d_A,    1, S2D_5pt, "double", OPS_READ),
  //         ops_arg_dat(d_Anew, 1, S2D_00,  "double", OPS_WRITE),
  //         ops_arg_reduce(h_err, 1, "double", OPS_MAX));

  
  
    //top bc
  //for (int i = 0; i < imax+2; i++)
  //  A[(jmax+1)*(imax+2)+i] = 0.0;

  // for (int j = 0; j < jmax+2; j++)
  // {
  //   A[(j)*(imax+2)+0] = sin(pi * j / (jmax+1));
  // }

  // for (int j = 0; j < imax+2; j++)
  // {
  //   A[(j)*(imax+2)+imax+1] = sin(pi * j / (jmax+1))*exp(-pi);
  // }

  // printf("Jacobi relaxation Calculation: %d x %d mesh\n", imax+2, jmax+2);

  // int iter = 0;

  // for (int i = 1; i < imax+2; i++)
  //   Anew[(0)*(imax+2)+i]   = 0.0;

  // for (int i = 1; i < imax+2; i++)
  //   Anew[(jmax+1)*(imax+2)+i] = 0.0;

  // for (int j = 1; j < jmax+2; j++)
  //   Anew[(j)*(imax+2)+0]   = sin(pi * j / (jmax+1));

  // for (int j = 1; j < jmax+2; j++)
  //   Anew[(j)*(imax+2)+jmax+1] = sin(pi * j / (jmax+1))*expf(-pi);

  int iter = 0;

    while ( error > tol && iter < iter_max )
  {
    ops_par_loop(apply_stencil, "apply_stencil", block, 2, interior_range, 
          ops_arg_dat(d_A,    1, S2D_5pt, "double", OPS_READ),
          ops_arg_dat(d_Anew, 1, S2D_00,  "double", OPS_WRITE),
          ops_arg_reduce(h_err, 1, "double", OPS_MAX));
    ops_reduction_result(h_err, &error);


    ops_par_loop(copy, "copy", block, 2, interior_range,
              ops_arg_dat(d_A,    1, S2D_00, "double", OPS_WRITE),
              ops_arg_dat(d_Anew, 1, S2D_00, "double", OPS_READ));
  
    if(iter % 10 == 0) printf("%5d, %0.6f\n", iter, error);        
    iter++;
  }

  // for (int i = 0; i < jmax; i++) {
  //   for (int j = 0; j < imax; j++) {
  //     fprintf(file, "%lf", A[j*imax+i]);
  //   }
  //   fprintf(file, "\n");
  // }
  ops_print_dat_to_txtfile(d_Anew, "results.txt");

  //fprintf(fptr, "hello, world\n");
  // while ( error > tol && iter < iter_max )
  // {
  //   error = 0.0;
  //   for( int j = 1; j < jmax+1; j++ )
  //   {
  //     for( int i = 1; i < imax+1; i++)
  //     {
  //       Anew[(j)*(imax+2)+i] = 0.25f * ( A[(j)*(imax+2)+i+1] + A[(j)*(imax+2)+i-1]
  //           + A[(j-1)*(imax+2)+i] + A[(j+1)*(imax+2)+i]);
  //       error = fmax( error, fabs((float)(Anew[(j)*(imax+2)+i]-A[(j)*(imax+2)+i])));
  //     }
  //   }

  //   for( int j = 1; j < jmax+1; j++ )
  //   {
  //     for( int i = 1; i < imax+1; i++)
  //     {
  //       A[(j)*(imax+2)+i] = Anew[(j)*(imax+2)+i];    
  //     }
  //   }
  //   if(iter % 10 == 0) printf("%5d, %0.6f\n", iter, (float)error);        
  //   iter++;
  // }

  ops_printf("%5d, %0.6f\n", iter, error);

  double err_diff = fabs((100.0*(error/2.421354960840227e-03))-100.0);
  printf("Total error is within %3.15E %% of the expected error\n",err_diff);
  if(err_diff < 0.001)
    printf("This run is considered PASSED\n");
  else
    printf("This test is considered FAILED\n");

  //fclose(file); // Close the file
  //fclose(fptr);
  //Finalising the OPS library
  ops_exit();

  free(A);
  free(Anew);

  
  return 0;
}

