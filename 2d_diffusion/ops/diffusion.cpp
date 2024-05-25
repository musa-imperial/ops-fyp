#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


double dt = 0.001;
double T = 1;
int    Nx = 40-2;
int    Ny = 40-2;
double mu = 0.1;

double a = 0.1;

double dx = 1.0;
double dy = 1.0;
double Lx = dx*(Nx+2-1);
double Ly = dy*(Ny+2-1);

// double Lx = 1;
// double Ly = 1;
// double dx = Lx / (Nx+2-1);
// double dy = Ly / (Ny+2-1);

double h = dx;

double hmudt=mu/h/h*dt;
double a0 = 1-4*hmudt;


#define OPS_2D
#include <ops_seq_v2.h>
#include "diffusion_kernels.h"


int main(int argc, const char** argv)
{
  ops_init(argc, argv,1);
	// Tracking progress to print to terminal
	//int prog = 0;
	
	// Precomputing end index of the grid array
	
  double *tmp = NULL;

  double ct0, ct1, et0, et1;

  double l2error = 0;
  double max_error = 0;
  
//   printf("T = %s \n", argv[1]);
//   double T = atof(argv[1]);

//   double T = 10;

  ops_block block = ops_decl_block(2, "my_grid");
  //ops_block v_block = ops_decl_block(2, "my_grid");

  int size[] = {Nx, Ny};
  int base[] = {0,0};
  int d_m[] =  {-1,-1};
  int d_p[] =  {1,1};

  ops_dat d_u    = ops_decl_dat(block, 1, size, base,
                               d_m, d_p, tmp,    "double", "u");

  ops_dat d_u_calc    = ops_decl_dat(block, 1, size, base,
                               d_m, d_p, tmp,    "double", "u_calc");

  ops_decl_const("dt",1,"double",&dt);
  ops_decl_const("T",1,"double",&T);
  ops_decl_const("Nx",1,"int",&Nx);
  ops_decl_const("Ny",1,"int",&Ny);
  ops_decl_const("mu",1,"double",&mu);

  ops_decl_const("a",1,"double",&a);

  ops_decl_const("dx",1,"double",&dx);
  ops_decl_const("dy",1,"double",&dy);
  ops_decl_const("Lx",1,"double",&Lx);
  ops_decl_const("Ly",1,"double",&Ly);
  ops_decl_const("h",1,"double",&h);
  ops_decl_const("hmudt",1,"double",&hmudt);
  ops_decl_const("a0",1,"double",&a0);

  //ops_printf("\nLx = %lf\n",Lx);


  int s2d_00[] = {0,0};
  ops_stencil S2D_00 = ops_decl_stencil(2,1,s2d_00,"0,0");

  int s2d_interior[] = {0,0, 1,0, -1,0, 0,1, 0,-1};
  ops_stencil S2D_INTERIOR = ops_decl_stencil(2,5,s2d_interior,"interior");

  //reduction_handle
  ops_reduction h_l2err = ops_decl_reduction_handle(sizeof(double), "double", "l2error");  

  
  ops_reduction h_maxerr = ops_decl_reduction_handle(sizeof(double), "double", "max_error"); 


  ops_partition("");


  
  //loop ranges

  int interior[] = {0, Nx, 0, Ny};

  int all[] = {-1, Nx+1, -1, Ny+1};

  ops_timers(&ct0, &et0);
  //int all[] = {0, Nx, 0, Ny};

  ops_par_loop(set_zero, "set_zero", block, 2, all,
      ops_arg_dat(d_u, 1, S2D_00, "double", OPS_WRITE));

  ops_par_loop(set_zero, "set_zero", block, 2, all,
    ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE));


  //set initial condition
  ops_par_loop(initial_condition, "initial_condition", block, 2, all,
        ops_arg_dat(d_u,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_idx());


  for (double t = 0.0; t < T; t += dt) {


    ops_par_loop(interior_stencil_u, "interior_stencil_u", block, 2, interior,
        ops_arg_dat(d_u,    1,   S2D_INTERIOR, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(copy, "copy", block, 2, all,
        ops_arg_dat(d_u,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_READ));
    

  }

  ops_timers(&ct1, &et1);
  ops_timing_output(std::cout);

    ops_par_loop(solution_check, "solution_check", block, 2, interior,
        ops_arg_dat(d_u,    1, S2D_00, "double", OPS_READ),
        ops_arg_idx(),
        ops_arg_reduce(h_l2err, 1, "double", OPS_INC),
        ops_arg_reduce(h_maxerr, 1, "double", OPS_MAX));

    ops_reduction_result(h_l2err, &l2error);
    ops_reduction_result(h_maxerr, &max_error);


  ops_print_dat_to_txtfile(d_u, "u_check.txt");

  //ops_printf("L2 Error: %0.6f\n", l2error); 
  //ops_printf("Max percentage Error: %0.6e\n", max_error*100); 

  ops_printf("\n%lf %le",et1-et0, max_error*100);

  //Finalising the OPS library
  
  ops_exit();
  free(tmp);
  return 0;
}

