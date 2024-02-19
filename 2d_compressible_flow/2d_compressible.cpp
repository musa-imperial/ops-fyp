#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>



#define OPS_2D
#include <ops_seq_v2.h>
#include "2d_compressible_kernels.h"
#include "utils.h"
#include "derivatives.h"

#include "global_params.h"
#include "global_ops_vars.h"
//#include "global_decl.h"

void build_datasets();
void initl();

int main(int argc, const char** argv)
{
  ops_init(argc, argv,1);

    //int i,j,itemp,k,n,nxm,iread,ni,nj,isave,longueur,imodulo;
    double ct0, ct1, et0, et1;

    build_datasets();

    // Global variables

    ops_decl_const("xlx", 1, "double", &xlx);
    ops_decl_const("yly", 1, "double", &yly);
    ops_decl_const("CFL", 1, "double", &CFL);
    ops_decl_const("dlx", 1, "double", &dlx);
    ops_decl_const("dly", 1, "double", &dly);
    ops_decl_const("dx", 1, "double", &dx);
    ops_decl_const("xmu", 1, "double", &xmu);
    ops_decl_const("xkt", 1, "double", &xkt);
    ops_decl_const("um0", 1, "double", &um0);
    ops_decl_const("vm0", 1, "double", &vm0);
    ops_decl_const("tm0", 1, "double", &tm0);

    ops_decl_const("xba", 1, "double", &xba);
    ops_decl_const("gma", 1, "double", &gma);
    ops_decl_const("chp", 1, "double", &chp);
    ops_decl_const("eta", 1, "double", &eta);
    ops_decl_const("uu0", 1, "double", &uu0);
    ops_decl_const("dlt", 1, "double", &dlt);
    ops_decl_const("um", 1, "double", &um);
    ops_decl_const("vm", 1, "double", &vm);
    ops_decl_const("tm", 1, "double", &tm);
    ops_decl_const("x", 1, "double", &x);
    ops_decl_const("y", 1, "double", &y);
    ops_decl_const("dy", 1, "double", &dy);

    //Square cylinder
    ops_decl_const("d", 1, "double", &d);
    ops_decl_const("radius", 1, "double", &radius);

    //initl
    
    ops_decl_const("pi", 1, "double", &pi);
    ops_decl_const("ct6", 1, "double", &ct6);


    ops_decl_const("roi", 1, "double", &roi);
    ops_decl_const("cci", 1, "double", &cci);
    ops_decl_const("tpi", 1, "double", &tpi);
    ops_decl_const("chv", 1, "double", &chv);

    ops_decl_const("udx", 1, "double", &udx);
    ops_decl_const("udy", 1, "double", &udy);

    ops_decl_const("udxx", 1, "double", &udxx);
    ops_decl_const("udyy", 1, "double", &udyy);

    ops_decl_const("utt", 1, "double", &utt);
    ops_decl_const("qtt", 1, "double", &qtt);
    


    initl();
    
    CFL = 0.25;
    dlt = CFL*dlx;  //changed from dx to dlx
    ops_update_const( "CFL", 1, "double", &CFL);
    ops_update_const( "dlt", 1, "double", &dlt);
    ops_printf("\nThe time step of the simulation is %lf\n",dlt);

    int all[] = {-1, nx+1, -1, ny+1};
    // //nxm = nx in original code
    ops_par_loop(average, "average", block, 2, all,
        ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_READ),
        ops_arg_reduce(h_um0, 1, "double", OPS_INC));
    ops_reduction_result(h_um0, &um0);

    // ops_par_loop(average, "average", block, 2, all,
    //     ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ),
    //     ops_arg_reduce(h_vm0, 1, "double", OPS_INC));
    // ops_reduction_result(h_vm0, &vm0);

    // ops_par_loop(average, "average", block, 2, all,
    //     ops_arg_dat(d_scp,    1, S2D_00, "double", OPS_READ),
    //     ops_arg_reduce(h_tm0, 1, "double", OPS_INC));
    // ops_reduction_result(h_tm0, &tm0);

    // um0 = um0/nx/ny;
    // vm0 = vm0/nx/ny;
    // tm0 = tm0/nx/ny;

    // ops_update_const( "um0", 1, "double", &um0);
    // ops_update_const( "vm0", 1, "double", &vm0);
    // ops_update_const( "tm0", 1, "double", &tm0);

    // ops_printf("\nThe average u velocity at t = 0 is: %lf\n",um0);
    // ops_printf("\nThe average v velocity at t = 0 is: %lf\n",vm0);
    // ops_printf("\nThe average T temperature at t = 0 is: %lf\n",tm0);


    //Beginning of the time loop

    for (int t=0; t<nt; t++) {
        //fluxx
        //derix(rou,nx,ny,tb1,xlx)
        // ops_par_loop(derix, "derix", block, 2, left,
        // ops_arg_dat(d_rou,    1, S2D_DX, "double", OPS_READ),
        // ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
        // ops_arg_gbl(&udx, 1, "double", OPS_READ));

        // ops_par_loop(derix, "derix", block, 2, right,
        // ops_arg_dat(d_rou,    1, S2D_DX, "double", OPS_READ),
        // ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
        // ops_arg_gbl(&udx, 1, "double", OPS_READ));

        // ops_par_loop(derix, "derix", block, 2, top,
        // ops_arg_dat(d_rou,    1, S2D_DX, "double", OPS_READ),
        // ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
        // ops_arg_gbl(&udx, 1, "double", OPS_READ));

        // ops_par_loop(derix, "derix", block, 2, bottom,
        // ops_arg_dat(d_rou,    1, S2D_DX, "double", OPS_READ),
        // ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
        // ops_arg_gbl(&udx, 1, "double", OPS_READ));

        // ops_par_loop(derix, "derix", block, 2, interior,
        // ops_arg_dat(d_rou,    1, S2D_DX, "double", OPS_READ),
        // ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
        // ops_arg_gbl(&udx, 1, "double", OPS_READ));



//         //deriy(rov,nx,ny,tb2,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_rov,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_rov,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_rov,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_rov,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_rov,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         //fluxx1
//         ops_par_loop(fluxx1, "fluxx1", block, 2, all,
//         ops_arg_dat(d_fro,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_READ));

//         //fluxx2
//         ops_par_loop(fluxx2, "fluxx2", block, 2, all,
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_rou,    1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_READ));
        
//         ops_par_loop(fluxx2, "fluxx2", block, 2, all,
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_rou,    1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ));

//         //derix(pre,nx,ny,tb3,xlx)
//         ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(d_pre,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(d_pre,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(d_pre,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(d_pre,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(d_pre,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         //derix(tb1,nx,ny,tb4,xlx)
//         ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         //deriy(tb2,nx,ny,tb5,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         //derxx(uuu,nx,ny,tb6,xlx)
//         ops_par_loop(derxx, "derxx", block, 2, left,
//         ops_arg_dat(d_uuu,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, right,
//         ops_arg_dat(d_uuu,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, top,
//         ops_arg_dat(d_uuu,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, bottom,
//         ops_arg_dat(d_uuu,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, interior,
//         ops_arg_dat(d_uuu,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         //deryy(uuu,nx,ny,tb7,yly)
//         ops_par_loop(deryy, "deryy", block, 2, left,
//         ops_arg_dat(d_uuu,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, right,
//         ops_arg_dat(d_uuu,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, top,
//         ops_arg_dat(d_uuu,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, bottom,
//         ops_arg_dat(d_uuu,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, interior,
//         ops_arg_dat(d_uuu,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         //derix(vvv,nx,ny,tb8,xlx)
//         ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(d_vvv,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(d_vvv,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(d_vvv,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(d_vvv,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(d_vvv,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         //deriy(tb8,nx,ny,tb9,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         //!fluxx3
//         ops_par_loop(fluxx3, 'fluxx3', block, 2, all,,
//         ops_arg_dat(d_tba, 1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_tb6, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb7, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb9, 1, S2D_00, "double", OPS_READ));

//         //!fluxx4
//         ops_par_loop(fluxx4, 'fluxx4', block, 2, all,,
//         ops_arg_dat(d_fru, 1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_tb3, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb4, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb5, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tba, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_eps, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_uuu, 1, S2D_00, "double", OPS_READ));

//         //fluxx2
//         ops_par_loop(fluxx2, "fluxx2", block, 2, all,
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_rou,    1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ));
        
//         ops_par_loop(fluxx2, "fluxx2", block, 2, all,
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_rov,    1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ));

// /////////////////////////////////////////////////
//         //deriy(pre,nx,ny,tb3,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_pre,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_pre,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_pre,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_pre,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_pre,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb3,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         //derix(tb1,nx,ny,tb4,xlx)
//         ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(d_tb1,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb4,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         //deriy(tb2,nx,ny,tb5,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_tb2,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb5,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         //derxx(vvv,nx,ny,tb6,xlx)
//         ops_par_loop(derxx, "derxx", block, 2, left,
//         ops_arg_dat(d_vvv,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, right,
//         ops_arg_dat(d_vvv,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, top,
//         ops_arg_dat(d_vvv,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, bottom,
//         ops_arg_dat(d_vvv,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         ops_par_loop(derxx, "derxx", block, 2, interior,
//         ops_arg_dat(d_vvv,    1, S2D_DXX, "double", OPS_READ),
//         ops_arg_dat(d_tb6,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udxx, 1, "double", OPS_READ));

//         //deryy(vvv,nx,ny,tb7,yly)
//         ops_par_loop(deryy, "deryy", block, 2, left,
//         ops_arg_dat(d_vvv,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, right,
//         ops_arg_dat(d_vvv,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, top,
//         ops_arg_dat(d_vvv,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, bottom,
//         ops_arg_dat(d_vvv,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         ops_par_loop(deryy, "deryy", block, 2, interior,
//         ops_arg_dat(d_vvv,    1, S2D_DYY, "double", OPS_READ),
//         ops_arg_dat(d_tb7,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udyy, 1, "double", OPS_READ));

//         //derix(uuu,nx,ny,tb8,xlx)
//         ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(d_uuu,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(d_uuu,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(d_uuu,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(d_uuu,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(d_uuu,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb8,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         //deriy(tb8,nx,ny,tb9,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_tb8,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb9,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         //!fluxx3
//         ops_par_loop(fluxx3, 'fluxx3', block, 2, all,,
//         ops_arg_dat(d_tbb, 1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_tb7, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb6, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb9, 1, S2D_00, "double", OPS_READ));

//         //!fluxx4
//         ops_par_loop(fluxx4, 'fluxx4', block, 2, all,,
//         ops_arg_dat(d_frv, 1, S2D_00, "double", OPS_WRITE),
//         ops_arg_dat(d_tb3, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb4, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tb5, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_tbb, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_eps, 1, S2D_00, "double", OPS_READ),
//         ops_arg_dat(d_vvv, 1, S2D_00, "double", OPS_READ));

//         //Equation for the temperature

//         //derix(scp,nx,ny,tb1,xlx)
//         ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(d_scp,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(d_scp,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(d_scp,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(d_scp,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(d_scp,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(d_tb1,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         //deriy(scp,nx,ny,tb2,yly)
//         ops_par_loop(deriy, "deriy", block, 2, left,
//         ops_arg_dat(d_scp,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, right,
//         ops_arg_dat(d_scp,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, top,
//         ops_arg_dat(d_scp,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, bottom,
//         ops_arg_dat(d_scp,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));

//         ops_par_loop(deriy, "deriy", block, 2, interior,
//         ops_arg_dat(d_scp,    1, S2D_DY, "double", OPS_READ),
//         ops_arg_dat(d_tb2,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udy, 1, "double", OPS_READ));




    }
    
    

    
    
    /*
  
    ***********************************************************
    ***********************************************************
    
    */


    ops_print_dat_to_txtfile(d_eps, "eps.txt");
    ops_print_dat_to_txtfile(d_vvv, "vvv.txt");
    ops_print_dat_to_txtfile(d_eee, "eee.txt");
    //ops_print_dat_to_txtfile(d_v, "v_check.txt");


    ops_timers(&ct1, &et1);
    ops_timing_output(std::cout);

    //ops_printf("\nTotal Wall time %lf\n",et1-et0);

    //Finalising the OPS library
  
    ops_exit();

    //free(temp);
 

    return 0;
}

