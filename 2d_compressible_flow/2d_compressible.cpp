#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


//double dt = 0.001;
//double T = 10;
//printf("\n-- Solving the problem up to time T = %.2f with a time-step (dt) of %.2f and spacial step (h) of %d --\n", T, dt, 1);
int nx = 129-2;
int ny = 129-2;
int nt = 10000;
int ns = 3;
int nf = 3;
int mx = nf*nx;
int my = nf*ny;  

//Square cylinder global variables
double d = 1.0;
double radius = d/2;

double xlx = 4*d;
double yly = 4*d;
double dlx = xlx/nx;
double dly = yly/ny;
double CFL, xmu, xkt, dx;
 
double um0 = 0;
double vm0 = 0; 
double tm0 = 0;

double xba, gma, chp, eta, uu0, dlt, um, vm, tm, x, y, dy;

//FROM INITL SUBROUTINE

double pi=acos(-1.);

//param 
double roi,cci,tpi,chv;

//initl params
double ct6;



#define OPS_2D
#include <ops_seq_v2.h>
#include "2d_compressible_kernels.h"
#include "initial_conditions.h"
#include "utils.h"


int main(int argc, const char** argv)
{
  ops_init(argc, argv,1);
	// Tracking progress to print to terminal
	//int prog = 0;
	
	// Precomputing end index of the grid array
  //printf("\n-- Solving the problem up to time T = %.2f with a time-step (dt) of %.2f and spacial step (h) of %d --\n", T, dt, 1);
	
//   double *u = NULL;
//   double *u_calc= NULL;
//   double *v= NULL;
//   double *v_calc= NULL;

    double *uuu = NULL;
    double *vvv = NULL;
    double *rho = NULL;
    double *eee = NULL;
    double *pre = NULL;
    double *tmp = NULL;
    double *rou = NULL;
    double *rov = NULL;
    double *wz = NULL;
    double *tuu = NULL;
    double *tvv = NULL;

    double *roe = NULL;
    double *tb1 = NULL;
    double *tb2 = NULL;
    double *tb3 = NULL;
    double *tb4 = NULL;
    double *tb5 = NULL;
    double *tb6 = NULL;
    double *tb7 = NULL;
    double *tb8 = NULL;
    double *tb9 = NULL;

    double *tba = NULL;
    double *tbb = NULL;
    double *fro = NULL;
    double *fru = NULL;
    double *frv = NULL;
    double *fre = NULL;
    double *gro = NULL;
    double *gru = NULL;
    double *grv = NULL;
    double *gre = NULL;
    double *rot = NULL;
    double *eps = NULL;
    double *ftp = NULL;
    double *gtp = NULL;
    double *scp = NULL;

    double *xx = NULL;
    double *yy = NULL;
    double *tf = NULL;

    double *coef = NULL;

    int i,j,itemp,k,n,nxm,iread,ni,nj,isave,longueur,imodulo;
    // double* xlx,yly,dlx,xmu,xkt,um0,vm0,tm0;
    // double* xba,gma,chp,eta,uu0,um,vm,tm,x,y,dy;

    // double d = 1.0;
    // double radius = d/2;

    double ct0, ct1, et0, et1;


    // Blocks

    ops_block block = ops_decl_block(2, "my_grid");

    // Block params
    int size[] = {nx, ny};
    int tf_size[] = {mx, my};
    int base[] = {0,0};
    int d_m[] =  {-1,-1};
    int d_p[] =  {1,1};


    // Datasets

    ops_dat d_uuu  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, uuu,  "double", "uuu");
    ops_dat d_vvv  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, vvv,  "double", "vvv");
    ops_dat d_rho  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, rho,  "double", "rho");
    ops_dat d_eee  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, eee,  "double", "eee");
    ops_dat d_pre  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, pre,  "double", "pre");
    ops_dat d_tmp  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tmp,  "double", "tmp");
    ops_dat d_rou  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, rou,  "double", "rou");
    ops_dat d_rov  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, rov,  "double", "rov");
    ops_dat d_wz   = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, wz,   "double", "wz");
    ops_dat d_tuu  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tuu,  "double", "tuu");
    ops_dat d_tvv  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tvv,  "double", "tvv");

    ops_dat d_roe  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, roe,  "double", "roe");
    ops_dat d_tb1  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb1,  "double", "tb1");
    ops_dat d_tb2  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb2,  "double", "tb2");
    ops_dat d_tb3  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb3,  "double", "tb3");
    ops_dat d_tb4  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb4,  "double", "tb4");
    ops_dat d_tb5  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb5,  "double", "tb5");
    ops_dat d_tb6  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb6,  "double", "tb6");
    ops_dat d_tb7  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb7,  "double", "tb7");
    ops_dat d_tb8  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb8,  "double", "tb8");
    ops_dat d_tb9  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tb9,  "double", "tb9");

    ops_dat d_tba  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tba,  "double", "tba");
    ops_dat d_tbb  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tbb,  "double", "tbb");
    ops_dat d_fro  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, fro,  "double", "fro");
    ops_dat d_fru  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, fru,  "double", "fru");
    ops_dat d_frv  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, frv,  "double", "frv");
    ops_dat d_fre  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, fre,  "double", "fre");
    ops_dat d_gro  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, gro,  "double", "gro");
    ops_dat d_gru  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, gru,  "double", "gru");
    ops_dat d_grv  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, grv,  "double", "grv");
    ops_dat d_gre  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, gre,  "double", "gre");
    ops_dat d_rot  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, rot,  "double", "rot");
    ops_dat d_eps  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, eps,  "double", "eps");
    ops_dat d_ftp  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, ftp,  "double", "ftp");
    ops_dat d_gtp  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, gtp,  "double", "gtp");
    ops_dat d_scp  = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, scp,  "double", "scp");

    ops_dat d_xx   = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, xx,   "double", "xx");
    ops_dat d_yy   = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, yy,   "double", "yy");
    ops_dat d_tf   = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, tf,   "double", "tf");

    ops_dat d_coef = ops_decl_dat(block, 1, size, base,
                                d_m, d_p, coef, "double", "coef");

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
    
    // Ranges

    int all[] = {-1, nx+1, -1, ny+1};


    // Stencils

    int s2d_00[] = {0,0};
    ops_stencil S2D_00 = ops_decl_stencil(2,1,s2d_00,"0,0");

    // Reduction handlers

    //Reduction handle
  ops_reduction h_um0 = ops_decl_reduction_handle(sizeof(double), "double", "um0");
  ops_reduction h_vm0 = ops_decl_reduction_handle(sizeof(double), "double", "vm0");

    // Initialisation of the variables

    //param(xlx,yly,xmu,xba,gma,chp,roi,cci,d,tpi,chv,uu0);
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
      // Initialisation end
    
    CFL = 0.25;
    dlt = CFL*dlx;  //changed from dx to dlx
    ops_update_const( "CFL", 1, "double", &CFL);
    ops_update_const( "dlt", 1, "double", &dlt);
    ops_printf("\nThe time step of the simulation is %lf\n",dlt);

    //nxm = nx in original code
    ops_par_loop(average, "average", block, 2, all,
        ops_arg_dat(d_uuu,    1, S2D_00, "double", OPS_READ),
        ops_arg_reduce(h_um0, 1, "double", OPS_INC));


    ops_par_loop(average, "average", block, 2, all,
        ops_arg_dat(d_vvv,    1, S2D_00, "double", OPS_READ),
        ops_arg_reduce(h_vm0, 1, "double", OPS_INC));

    // ops_par_loop(average, "average", block, 2, all,
    //     ops_arg_dat(d_scp,    1, S2D_00, "double", OPS_READ),
    //     ops_arg_reduce(tm0, 1, "double", OPS_INC));

    um0 = um0/nx/ny;
    vm0 = vm0/nx/ny;
    tm0 = tm0/nx/ny;

    ops_update_const( "um0", 1, "double", &um0);
    ops_update_const( "vm0", 1, "double", &vm0);
    ops_update_const( "tm0", 1, "double", &tm0);

    ops_printf("\nThe average v velocity at t = 0 is: %lf\n",vm0);
    
    

    
    
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

    ops_printf("\nTotal Wall time %lf\n",et1-et0);

    //Finalising the OPS library
  
    ops_exit();

    free(uuu);
    free(vvv);
    free(rho);
    free(eee);
    free(pre);
    free(tmp);
    free(rou);
    free(rov);
    free(wz);
    free(tuu);
    free(tvv);

    free(roe);
    free(tb1);
    free(tb2);
    free(tb3);
    free(tb4);
    free(tb5);
    free(tb6);
    free(tb7);
    free(tb8);
    free(tb9);

    free(tba);
    free(tbb);
    free(fro);
    free(fru);
    free(frv);
    free(fre);
    free(gro);
    free(gru);
    free(grv);
    free(gre);
    free(rot);
    free(eps);
    free(ftp);
    free(gtp);
    free(scp);

    free(xx);
    free(yy);
    free(tf);

    return 0;
}

