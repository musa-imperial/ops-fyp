#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// OPS header file
#define OPS_2D
#include "ops_seq_v2.h"

#include "global_ops_vars.h"

#include "global_params.h"

void build_datasets() 
{
    // block
    block = ops_decl_block(2, "2D_grid");
    
    // Block params
    int size[] = {nx, ny};
    int tf_size[] = {mx, my};
    int base[] = {0,0};
    int d_m[] =  {-1,-1};
    int d_p[] =  {1,1};
    double* temp = NULL;

    //datasets
    d_uuu  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "uuu");
    d_vvv  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "vvv");
    d_rho  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "rho");
    d_eee  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "eee");
    d_pre  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "pre");
    d_tmp  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tmp");
    d_rou  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "rou");
    d_rov  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "rov");
    d_wz   = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "wz");
    d_tuu  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tuu");
    d_tvv  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tvv");

    d_roe  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "roe");
    d_tb1  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb1");
    d_tb2  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb2");
    d_tb3  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb3");
    d_tb4  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb4");
    d_tb5  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb5");
    d_tb6  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb6");
    d_tb7  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb7");
    d_tb8  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb8");
    d_tb9  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tb9");

    d_tba  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tba");
    d_tbb  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tbb");
    d_fro  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "fro");
    d_fru  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "fru");
    d_frv  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "frv");
    d_fre  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "fre");
    d_gro  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "gro");
    d_gru  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "gru");
    d_grv  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "grv");
    d_gre  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "gre");
    d_rot  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "rot");
    d_eps  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "eps");
    d_ftp  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "ftp");
    d_gtp  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "gtp");
    d_scp  = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "scp");

    d_xx   = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "xx");
    d_yy   = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "yy");
    d_tf   = ops_decl_dat(block, 1, size, base, d_m, d_p, temp,  "double", "tf");

    d_coef = ops_decl_dat(block, 1, tf_size, base, d_m, d_p, temp, "double", "coef");

    // Stencils

    int s2d_00[] = {0,0};
    int s2d_dx[] = {1,0, 0,-1};
    int s2d_dy[] = {0,1, -1,0};
    int s2d_dxx[] = {0,0, 1,0, 0,-1};
    int s2d_dyy[] = {0,0, 0,1, -1,0};

    S2D_00 = ops_decl_stencil(2,1,s2d_00,"0,0");

    S2D_DX = ops_decl_stencil(2,2,s2d_dx,"dx");

    S2D_DY = ops_decl_stencil(2,2,s2d_dy,"dy");

    S2D_DXX = ops_decl_stencil(2,3,s2d_dxx,"dxx");

    S2D_DYY = ops_decl_stencil(2,3,s2d_dyy,"dyy");

    // Reduction handlers

    h_um0 = ops_decl_reduction_handle(sizeof(double), "double", "um0");
    h_vm0 = ops_decl_reduction_handle(sizeof(double), "double", "vm0");
    h_tm0 = ops_decl_reduction_handle(sizeof(double), "double", "tm0");
}