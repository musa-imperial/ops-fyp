#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


double dt = 0.001;
double T = 100;
//printf("\n-- Solving the problem up to time T = %.2f with a time-step (dt) of %.2f and spacial step (h) of %d --\n", T, dt, 1);
int    Nx = 501;
int    Ny = 501;
double a = 0.75;
double b = 0.06;
double mu1 = 5.0;
double mu2 = 0.0;
double eps = 13.0;

double dx = 1.0;
double dy = 1.0;
double Lx = dx*(Nx-1);
double Ly = dy*(Ny-1);
double h = dx;

double hmu1dt=mu1/h/h*dt;
double hmu2dt=mu2/h/h*dt;
double div_a = 1/a;


#define OPS_2D
#include <ops_seq_v2.h>
#include "reaction_diffusion_kernels.h"


int main(int argc, const char** argv)
{
  ops_init(argc, argv,1);
	// Tracking progress to print to terminal
	//int prog = 0;
	
	// Precomputing end index of the grid array
  //printf("\n-- Solving the problem up to time T = %.2f with a time-step (dt) of %.2f and spacial step (h) of %d --\n", T, dt, 1);
	
  double *u = NULL;
  double *u_calc= NULL;
  double *v= NULL;
  double *v_calc= NULL;

  double ct0, ct1, et0, et1;
  
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
                               d_m, d_p, u,    "double", "u");

  ops_dat d_u_calc    = ops_decl_dat(block, 1, size, base,
                               d_m, d_p, u_calc,    "double", "u_calc");

  ops_dat d_v    = ops_decl_dat(block, 1, size, base,
                               d_m, d_p, v,    "double", "v");

  ops_dat d_v_calc    = ops_decl_dat(block, 1, size, base,
                               d_m, d_p, v_calc,    "double", "v_calc");

  ops_decl_const("dt",1,"double",&dt);
  ops_decl_const("T",1,"double",&T);
  ops_decl_const("Nx",1,"int",&Nx);
  ops_decl_const("Ny",1,"int",&Ny);
  ops_decl_const("a",1,"double",&a);
  ops_decl_const("b",1,"double",&b);
  ops_decl_const("mu1",1,"double",&mu1);
  ops_decl_const("mu2",1,"double",&mu2);
  ops_decl_const("eps",1,"double",&eps);

  ops_decl_const("dx",1,"double",&dx);
  ops_decl_const("dy",1,"double",&dy);
  ops_decl_const("Lx",1,"double",&Lx);
  ops_decl_const("Ly",1,"double",&Ly);
  ops_decl_const("h",1,"double",&h);

  ops_decl_const("hmu1dt",1,"double",&hmu1dt);
  ops_decl_const("hmu2dt",1,"double",&hmu2dt);
  ops_decl_const("div_a",1,"double",&div_a);
  //ops_decl_const("pi",1,"double",&pi);


  int s2d_00[] = {0,0};
  ops_stencil S2D_00 = ops_decl_stencil(2,1,s2d_00,"0,0");

  //corner bc stencils
  //point x = 0,  y = 0 stencil (bottom_left)
  int s2d_bottom_left[] = {0,0,   1,0,  0,1};
  ops_stencil S2D_BOTTOM_LEFT = ops_decl_stencil(2,3,s2d_bottom_left,"bottom_left");

  //point x = 0,  y = Ny-1 stencil (top_left)
  int s2d_top_left[] = {0,0,  1,0, 0,-1};
  ops_stencil S2D_TOP_LEFT = ops_decl_stencil(2,3,s2d_top_left,"top_left");

  //point x = Nx-1, y = 0 stencil (bottom_right)
  int s2d_bottom_right[] = {0,0,  -1,0, 0,1};
  ops_stencil S2D_BOTTOM_RIGHT = ops_decl_stencil(2,3,s2d_bottom_right,"bottom_right");

  //point x = Nx-1, y = Ny-1 stencil (top_right)
  int s2d_top_right[] = {0,0, -1,0,  0,-1};
  ops_stencil S2D_TOP_RIGHT = ops_decl_stencil(2,3,s2d_top_right,"top_right");

  //bc stencils
  int s2d_left[] = {0,0, 1,0, 0,1, 0,-1};
  ops_stencil S2D_LEFT = ops_decl_stencil(2,4,s2d_left,"left");

  int s2d_right[] = {0,0, -1,0, 0,1, 0,-1};
  ops_stencil S2D_RIGHT = ops_decl_stencil(2,4,s2d_right,"right");

  int s2d_top[] = {0,0, 1,0, -1,0, 0,-1};
  ops_stencil S2D_TOP = ops_decl_stencil(2,4,s2d_top,"top");

  int s2d_bottom[] = {0,0, 1,0, -1,0, 0,1};
  ops_stencil S2D_BOTTOM = ops_decl_stencil(2,4,s2d_bottom,"bottom");

  int s2d_interior[] = {0,0, 1,0, -1,0, 0,1, 0,-1};
  ops_stencil S2D_INTERIOR = ops_decl_stencil(2,5,s2d_interior,"interior");


  ops_partition("");


  
  //loop ranges
  int bottom_left[] = {-1, 0, -1, 0};

  int bottom_right[] = {Nx, Nx+1, -1, 0};

  int top_left[] = {-1, 0, Ny, Ny+1};

  int top_right[] = {Nx, Nx+1, Ny, Ny+1};

  int bottom[] = {0, Nx, -1, 0};

  int top[] = {0, Nx, Ny, Ny+1};

  int left[] = {-1, 0, 0, Ny};

  int right[] = {Nx, Nx+1, 0, Ny};

  int interior[] = {0, Nx, 0, Ny};

  int all[] = {-1, Nx+1, -1, Ny+1};

  ops_timers(&ct0, &et0);
  //int all[] = {0, Nx, 0, Ny};

//   ops_par_loop(set_zero, "set_zero", block, 2, bottom,
//       ops_arg_dat(d_u, 1, S2D_00, "double", OPS_WRITE));

//   ops_par_loop(set_zero, "set_zero", block, 2, all,
//     ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE));

//   ops_par_loop(set_zero, "set_zero", block, 2, all,
//       ops_arg_dat(d_v, 1, S2D_00, "double", OPS_WRITE));

//   ops_par_loop(set_zero, "set_zero", block, 2, all,
//     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));


  //set initial condition
  ops_par_loop(u_initcond_stencil, "u_initcond_stencil", block, 2, all,
        ops_arg_dat(d_u,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_idx());

  ops_par_loop(v_initcond_stencil, "v_initcond_stencil", block, 2, all,
        ops_arg_dat(d_v,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_idx());


  for (double t = 0.0; t < T; t += dt) {

    //  Corner boundary conditions
    ops_par_loop(bottomleft_u, "bottomleft_u", block, 2, bottom_left,
        ops_arg_dat(d_u,    1,   S2D_BOTTOM_LEFT, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_BOTTOM_LEFT, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(bottomleft_v, "bottomleft_v", block, 2, bottom_left,
    //     ops_arg_dat(d_v,    1,   S2D_BOTTOM_LEFT, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(topleft_u, "topleft_u", block, 2, top_left,
        ops_arg_dat(d_u,    1,   S2D_TOP_LEFT, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_TOP_LEFT, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(topleft_v, "topleft_v", block, 2, top_left,
    //     ops_arg_dat(d_v,    1,   S2D_TOP_LEFT, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottomright_u, "bottomright_u", block, 2, bottom_right,
        ops_arg_dat(d_u,    1,   S2D_BOTTOM_RIGHT, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_BOTTOM_RIGHT, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(bottomright_v, "bottomright_v", block, 2, bottom_right,
    //     ops_arg_dat(d_v,    1,   S2D_BOTTOM_RIGHT, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(topright_u, "topright_u", block, 2, top_right,
        ops_arg_dat(d_u,    1,   S2D_TOP_RIGHT, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_TOP_RIGHT, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(topright_v, "topright_v", block, 2, top_right,
    //     ops_arg_dat(d_v,    1,   S2D_TOP_RIGHT, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    //

    ops_par_loop(left_bndcon_u, "left_bndcon_u", block, 2, left,
        ops_arg_dat(d_u,    1,   S2D_LEFT, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_LEFT, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(left_bndcon_v, "left_bndcon_v", block, 2, left,
    //     ops_arg_dat(d_v,    1,   S2D_LEFT, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(right_bndcon_u, "right_bndcon_u", block, 2, right,
        ops_arg_dat(d_u,    1,   S2D_RIGHT, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_RIGHT, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(right_bndcon_v, "right_bndcon_v", block, 2, right,
    //     ops_arg_dat(d_v,    1,   S2D_RIGHT, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ///
    ops_par_loop(top_bndcon_u, "top_bndcon_u", block, 2, top,
        ops_arg_dat(d_u,    1,   S2D_TOP, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_TOP, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(top_bndcon_v, "top_bndcon_v", block, 2, top,
    //     ops_arg_dat(d_v,    1,   S2D_TOP, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottom_bndcon_u, "bottom_bndcon_u", block, 2, bottom,
        ops_arg_dat(d_u,    1,   S2D_BOTTOM, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_BOTTOM, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(bottom_bndcon_v, "bottom_bndcon_v", block, 2, bottom,
    //     ops_arg_dat(d_v,    1,   S2D_BOTTOM, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(interior_stencil_u, "interior_stencil_u", block, 2, interior,
        ops_arg_dat(d_u,    1,   S2D_INTERIOR, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_INTERIOR, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));

    // ops_par_loop(interior_stencil_v, "interior_stencil_v", block, 2, interior,
    //     ops_arg_dat(d_v,    1,   S2D_INTERIOR, "double", OPS_READ),
    //     ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
    //     ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    // std::swap<double*>(u, u_calc);
	// std::swap<double*>(v, v_calc);
    ops_par_loop(copy, "copy", block, 2, all,
        ops_arg_dat(d_u,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_READ));
    
    ops_par_loop(copy, "copy", block, 2, all,
        ops_arg_dat(d_v,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_READ));
    

  }
  ops_print_dat_to_txtfile(d_u, "u_check.txt");
  ops_print_dat_to_txtfile(d_v, "v_check.txt");

  
  ops_timers(&ct1, &et1);
  ops_timing_output(std::cout);

  ops_printf("\nTotal Wall time %lf\n",et1-et0);

  //Finalising the OPS library

  //ops_printf("%lf", u[0]);
  
  ops_exit();
  free(u);
  free(u_calc);
  free(v);
  free(v_calc);
  return 0;
}

