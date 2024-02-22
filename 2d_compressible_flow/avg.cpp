#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>



#define OPS_2D
#include <ops_seq_v2.h>

#include "global_ops_vars.h"

#include "global_params.h"

#include "avg_kernel.h"
void avg_t0() {

    int all[] = {-1, nx+1, -1, ny+1};

    // ops_par_loop(average, "average", block, 2, all,
    //     ops_arg_dat(d_u1,    1, S2D_00, "double", OPS_READ),
    //     ops_arg_reduce(h_u1, 1, "double", OPS_INC));

    ops_par_loop(average, "average", block, 2, all,
        ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_READ),
        ops_arg_reduce(h_um0, 1, "double", OPS_INC));
    
    ops_par_loop(average, "average", block, 2, all,
        ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ),
        ops_arg_reduce(h_vm0, 1, "double", OPS_INC));
    ops_reduction_result(h_vm0, &vm0);

    ops_par_loop(average, "average", block, 2, all,
        ops_arg_dat(d_scp,    1, S2D_00, "double", OPS_READ),
        ops_arg_reduce(h_tm0, 1, "double", OPS_INC));
    ops_reduction_result(h_tm0, &tm0);

    um0 = um0/nx/ny;
    vm0 = vm0/nx/ny;
    tm0 = tm0/nx/ny;

    ops_update_const( "um0", 1, "double", &um0);
    ops_update_const( "vm0", 1, "double", &vm0);
    ops_update_const( "tm0", 1, "double", &tm0);

    ops_printf("\nThe average u velocity at t = 0 is: %lf\n",um0);
    ops_printf("\nThe average v velocity at t = 0 is: %lf\n",vm0);
    ops_printf("\nThe average T temperature at t = 0 is: %lf\n",tm0);

}