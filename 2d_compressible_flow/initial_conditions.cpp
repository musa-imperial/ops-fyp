
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// OPS header file
#define OPS_2D
#include "ops_seq_v2.h"

#include "global_ops_vars.h"

#include "initial_conditions_kernels.h"

void initl() {
    int bottom_left[] = {-1, 0, -1, 0};

    int bottom_right[] = {nx, nx+1, -1, 0};

    int top_left[] = {-1, 0, ny, ny+1};

    int top_right[] = {nx, nx+1, ny, ny+1};

    int bottom[] = {0, nx, -1, 0};

    int top[] = {0, nx, ny, ny+1};

    int left[] = {-1, 0, 0, ny};

    int right[] = {nx, nx+1, 0, ny};

    int interior[] = {0, nx, 0, ny};

    int all[] = {-1, nx+1, -1, ny+1};

    ops_par_loop(square_cylinder, "square_cylinder", block, 2, all,
        ops_arg_dat(d_eps,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_idx());

    ops_par_loop(init_condition_uuu, "init_condition_uuu", block, 2, all,
        ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(init_condition_vvv, "init_condition_vvv", block, 2, all,
        ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_idx());
    
    ops_par_loop(init_condition_tmp, "init_condition_tmp", block, 2, all,
        ops_arg_dat(d_tmp,    1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(init_condition_eee, "init_condition_eee", block, 2, all,
        ops_arg_dat(d_eee,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_tmp,    1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ));

    ops_par_loop(init_condition_rho, "init_condition_rho", block, 2, all,
        ops_arg_dat(d_rho,    1, S2D_00, "double", OPS_WRITE));

    ops_par_loop(init_condition_pre, "init_condition_pre", block, 2, all,
        ops_arg_dat(d_pre,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_rho,    1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_tmp,    1, S2D_00, "double", OPS_READ));

    ops_par_loop(init_condition_rou, "init_condition_rou", block, 2, all,
        ops_arg_dat(d_rou,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_rho,    1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_READ));

    ops_par_loop(init_condition_rov, "init_condition_rov", block, 2, all,
        ops_arg_dat(d_rov,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_rho,    1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ));

    ops_par_loop(init_condition_roe, "init_condition_roe", block, 2, all,
        ops_arg_dat(d_roe,    1, S2D_00, "double", OPS_WRITE),
        ops_arg_dat(d_rho,    1, S2D_00, "double", OPS_READ),
        ops_arg_dat(d_eee,    1, S2D_00, "double", OPS_READ));

    ops_par_loop(init_condition_scp, "init_condition_scp", block, 2, all,
        ops_arg_dat(d_scp,    1, S2D_00, "double", OPS_WRITE));

}