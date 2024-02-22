
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// OPS header file
#define OPS_2D
#include "ops_seq_v2.h"

#include "global_params.h"
#include "global_ops_vars.h"
#include "initial_conditions_kernel.h"


void initl() {
    double ren, pdl, mach;
    //FROM PARAMS SUBROUTINE

    ren=200.0; //REYNOLDS NUMBER
    mach=0.2; //MACH NUMBER
    pdl=0.7;  //PRANDTL NUMBER
    roi=1.0;
    cci=1.0;   //SOUND SPEED
    d=1.0;
    chp=1.0;
    gma=1.4;
        
    chv=chp/gma;
    xlx=4.0*d;     //DOMAIN SIZE X DIRECTION
    yly=4.0*d;     //DOMAIN SIZE Y DIRECTION
    uu0=mach*cci;
    xmu=roi*uu0*d/ren;
    xba=xmu*chp/pdl;
    tpi=cci*cci/(chp*(gma-1));

    //FROM INITL SUBROUTINE
    
    double epsi=0.1;
    //dlx=xlx/nx
    //dly=yly/ny
    double ct3=log(2.);
    double ct4=yly/2.;
    double ct5=xlx/2.;
    double ct6=(gma-1.)/gma;
    ops_update_const( "ct6", 1, "double", &ct6);
    double y=-ct4;
    double x=0.;
    double eta=0.1;
    eta=eta/2.;
    double radius=d/2.;
    double xkt=xba/(chp*roi);
    //double pi=acos(-1.);

    ops_update_const( "xlx", 1, "double", &xlx);
    ops_update_const( "yly", 1, "double", &yly);
    ops_update_const( "roi", 1, "double", &cci);
    ops_update_const( "d", 1, "double", &d);
    ops_update_const( "chp", 1, "double", &chp);
    ops_update_const( "gma", 1, "double", &gma);
    ops_update_const( "uu0", 1, "double", &uu0);
    ops_update_const( "xmu", 1, "double", &xmu);
    ops_update_const( "xba", 1, "double", &xba);
    ops_update_const( "tpi", 1, "double", &tpi);

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