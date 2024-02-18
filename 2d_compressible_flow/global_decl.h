//
//blocks
//

ops_block block;

//
//extern ops dats
//
ops_dat d_uuu;
ops_dat d_vvv;
ops_dat d_rho;
ops_dat d_eee;
ops_dat d_pre;
ops_dat d_tmp;
ops_dat d_rou;
ops_dat d_rov;
ops_dat d_wz;
ops_dat d_tuu;
ops_dat d_tvv;

ops_dat d_roe;
ops_dat d_tb1;
ops_dat d_tb2;
ops_dat d_tb3;
ops_dat d_tb4;
ops_dat d_tb5;
ops_dat d_tb6;
ops_dat d_tb7;
ops_dat d_tb8;
ops_dat d_tb9;

ops_dat d_tba;
ops_dat d_tbb;
ops_dat d_fro;
ops_dat d_fru;
ops_dat d_frv;
ops_dat d_fre;
ops_dat d_gro;
ops_dat d_gru;
ops_dat d_grv;
ops_dat d_gre;
ops_dat d_rot;
ops_dat d_eps;
ops_dat d_ftp;
ops_dat d_gtp;
ops_dat d_scp;

ops_dat d_xx;
ops_dat d_yy;
ops_dat d_tf;

ops_dat d_coef;


//
// Stencils
//
ops_stencil S2D_00;

ops_stencil S2D_DX;

ops_stencil S2D_DY;

ops_stencil S2D_DXX;

ops_stencil S2D_DYY;

//
// Reductions
//
ops_reduction h_um0;
ops_reduction h_vm0;
ops_reduction h_tm0;

//Global params

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

//derix, deriy, derxx, deryy params
double udx = 1/(dlx+dlx);
double udy = 1/(dly+dly);

double udxx = 1/(dlx*dlx);
double udyy = 1/(dly*dly);

//fluxx params
double utt=1/3;
double qtt=4/3;