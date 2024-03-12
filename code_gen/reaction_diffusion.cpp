
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

double dt = 0.001;
double T = 100;
int Nx = 251;
int Ny = 251;
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
double hmu1dt = mu1/h/h*dt;
double hmu2dt = mu2/h/h*dt;
double div_a = 1/a;

#define OPS_2D
#include <ops_seq_v2.h>
#include "reaction_diffusion_kernels.h" 
int main(int argc, const char** argv)
{
    ops_init(argc, argv,1);

    // block
    ops_block block = ops_decl_block(2, "2D_grid");
    
    // Block params
    int size[] = {Nx, Ny};
    int base[] = {0,0};
    int d_m[] =  {-1,-1};
    int d_p[] =  {1,1};
    double* temp = NULL;

    ops_dat d_u  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp, "double", "u");
    ops_dat d_v  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp, "double", "v");
    ops_dat d_u_calc  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp, "double", "u_calc");
    ops_dat d_v_calc  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp, "double", "v_calc");
    
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

    
    ops_par_loop(set_zero, "set_zero", block, 2, all,
            ops_arg_dat(d_u, 1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(set_zero, "set_zero", block, 2, all,
            ops_arg_dat(d_v, 1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(set_zero, "set_zero", block, 2, all,
            ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(set_zero, "set_zero", block, 2, all,
            ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE));


    ops_par_loop(u_initcond_stencil, "u_initcond_stencil", block, 2, all,
            ops_arg_dat(d_u,    1, S2D_00, "double", OPS_WRITE),
            ops_arg_idx());

    ops_par_loop(v_initcond_stencil, "v_initcond_stencil", block, 2, all,
            ops_arg_dat(d_v,    1, S2D_00, "double", OPS_WRITE),
            ops_arg_idx());

for (double t = 0; t < T; t += dt) {
    
    ops_par_loop(bottom_left_u, "bottom_left_u", block, 2, bottom_left,
                ops_arg_dat(d_u,    1,   S2D_BOTTOM_LEFT, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottom_left_v, "bottom_left_v", block, 2, bottom_left,
                ops_arg_dat(d_v,    1,   S2D_BOTTOM_LEFT, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottom_right_u, "bottom_right_u", block, 2, bottom_right,
                ops_arg_dat(d_u,    1,   S2D_BOTTOM_RIGHT, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottom_right_v, "bottom_right_v", block, 2, bottom_right,
                ops_arg_dat(d_v,    1,   S2D_BOTTOM_RIGHT, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(top_left_u, "top_left_u", block, 2, top_left,
                ops_arg_dat(d_u,    1,   S2D_TOP_LEFT, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(top_left_v, "top_left_v", block, 2, top_left,
                ops_arg_dat(d_v,    1,   S2D_TOP_LEFT, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(top_right_u, "top_right_u", block, 2, top_right,
                ops_arg_dat(d_u,    1,   S2D_TOP_RIGHT, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(top_right_v, "top_right_v", block, 2, top_right,
                ops_arg_dat(d_v,    1,   S2D_TOP_RIGHT, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(left_u, "left_u", block, 2, left,
                ops_arg_dat(d_u,    1,   S2D_LEFT, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(left_v, "left_v", block, 2, left,
                ops_arg_dat(d_v,    1,   S2D_LEFT, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(right_u, "right_u", block, 2, right,
                ops_arg_dat(d_u,    1,   S2D_RIGHT, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(right_v, "right_v", block, 2, right,
                ops_arg_dat(d_v,    1,   S2D_RIGHT, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(top_u, "top_u", block, 2, top,
                ops_arg_dat(d_u,    1,   S2D_TOP, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(top_v, "top_v", block, 2, top,
                ops_arg_dat(d_v,    1,   S2D_TOP, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottom_u, "bottom_u", block, 2, bottom,
                ops_arg_dat(d_u,    1,   S2D_BOTTOM, "double", OPS_READ),
                ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(bottom_v, "bottom_v", block, 2, bottom,
                ops_arg_dat(d_v,    1,   S2D_BOTTOM, "double", OPS_READ),
                ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
                ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    
    ops_par_loop(interior_stencil_u, "interior_stencil_u", block, 2, interior,
        ops_arg_dat(d_u,    1,   S2D_INTERIOR, "double", OPS_READ),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v,    1,   S2D_00, "double", OPS_READ));

    ops_par_loop(interior_stencil_v, "interior_stencil_v", block, 2, interior,
        ops_arg_dat(d_v,    1,   S2D_INTERIOR, "double", OPS_READ),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_u,    1,   S2D_00, "double", OPS_READ));

    
    
    ops_par_loop(copy, "copy", block, 2, all,
        ops_arg_dat(d_u,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_u_calc, 1, S2D_00, "double", OPS_READ));

    ops_par_loop(copy, "copy", block, 2, all,
        ops_arg_dat(d_v,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_v_calc, 1, S2D_00, "double", OPS_READ));

    }
    
        ops_print_dat_to_txtfile(d_u, "u_results.txt");   
        
        ops_print_dat_to_txtfile(d_v, "v_results.txt");   
        
    //Finalising the OPS library

    ops_exit();
    free(temp);
    
    return 0;
}
