// void param(double *xlx, double *yly, double *xmu, double *xba, double *gma,
//                   double *chp, double *roi, double *cci, double *d, double *tpi,
//                   double *chv, double *uu0) {
//     double ren, pdl, mach;

//     ren=200.0; //REYNOLDS NUMBER
//     mach=0.2; //MACH NUMBER
//     pdl=0.7;  //PRANDTL NUMBER
//     roi=1.0;
//     cci=1.0;   //SOUND SPEED
//     d=1.0;
//     chp=1.0;
//     gma=1.4;
        
//     chv=chp/gma;
//     xlx=4.0*d;     //DOMAIN SIZE X DIRECTION
//     yly=4.0*d;     //DOMAIN SIZE Y DIRECTION
//     uu0=mach*cci;
//     xmu=roi*uu0*d/ren;
//     xba=xmu*chp/pdl;
//     tpi=cci**2/(chp*(gma-1));

// }

// void derix(ops_dat phi, ops_dat dfi) {
//     ops_par_loop(derix, "derix", block, 2, left,
//         ops_arg_dat(phi,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(dfi,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, right,
//         ops_arg_dat(phi,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(dfi,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, top,
//         ops_arg_dat(phi,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(dfi,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, bottom,
//         ops_arg_dat(phi,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(dfi,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));

//         ops_par_loop(derix, "derix", block, 2, interior,
//         ops_arg_dat(phi,    1, S2D_DX, "double", OPS_READ),
//         ops_arg_dat(dfi,    1, S2D_00, "double", OPS_WRITE),
//         ops_arg_gbl(&udx, 1, "double", OPS_READ));
// }