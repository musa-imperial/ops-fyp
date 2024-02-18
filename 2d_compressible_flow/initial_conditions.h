
#include "global_decl.h"

void set_zero(ACC<double> &A) {
    A(0,0) = 0.0;
}

void copy(ACC<double> &A, const ACC<double> &Anew) {
    A(0,0) = Anew(0,0);
}

// void initl(ACC<double> &A, const int *idx) {

//     if (idx[1] > Ly/2) {
//         A(0, 0) = 1;
//     }
//     else {
//         A(0, 0) = 0;
//     }
    
// }

void square_cylinder(ACC<double> &eps, const int *idx) {
    if ((idx[0]*dlx-xlx/2)*(idx[0]*dlx-xlx/2) + (idx[1]*dly-yly/2)*(idx[1]*dly-yly/2) < radius*radius) {
        eps(0,0) = 1;
    }
    else {
        eps(0,0) = 0;
    }
}

void init_condition_uuu(ACC<double> &uuu) {
    uuu(0,0) = uu0;
}

void init_condition_vvv(ACC<double> &vvv, const int *idx) {
    vvv(0,0) = 0.01*(sin(4.*pi*idx[0]*dlx/xlx)
                +sin(7.*pi*idx[0]*dlx/xlx))
                *exp(-(idx[1]*dly-yly/2.)*(idx[1]*dly-yly/2.));
}

void init_condition_tmp(ACC<double> &tmp) {
    tmp(0,0) = tpi;
}

void init_condition_eee(ACC<double> &eee, const ACC<double> &tmp, const ACC<double> &uuu, const ACC<double> &vvv) {
    eee(0,0)=chv*tmp(0,0)+0.5*(uuu(0, 0)*uuu(0,0)+vvv(0,0)*vvv(0,0));
}

void init_condition_rho(ACC<double> &rho) {
    rho(0,0)=roi;
}

void init_condition_pre(ACC<double> &pre, const ACC<double> &rho, const ACC<double> &tmp) {
    pre(0,0)=rho(0,0)*ct6*chp*tmp(0,0);
}

void init_condition_rou(ACC<double> &rou, const ACC<double> &rho, const ACC<double> &uuu) {
    rou(0,0)=rho(0,0)*uuu(0,0);
}

void init_condition_rov(ACC<double> &rov, const ACC<double> &rho, const ACC<double> &vvv) {
    rov(0,0)=rho(0,0)*vvv(0,0);
}

void init_condition_roe(ACC<double> &roe, const ACC<double> &rho, const ACC<double> &eee) {
    roe(0,0)=rho(0,0)*eee(0,0);
}

void init_condition_scp(ACC<double> &scp) {
    scp(0,0)=1;
}
 

//   do j=1,ny
//      do i=1,nx
//         uuu(i,j)=uu0
//         vvv(i,j)=0.01*(sin(4.*pi*i*dlx/xlx)&
//              +sin(7.*pi*i*dlx/xlx))*&
//              exp(-(j*dly-yly/2.)**2)
//         tmp(i,j)=tpi
//         eee(i,j)=chv*tmp(i,j)+0.5*(uuu(i,j)*uuu(i,j)+vvv(i,j)*vvv(i,j))
//         rho(i,j)=roi
//         pre(i,j)=rho(i,j)*ct6*chp*tmp(i,j)
//         rou(i,j)=rho(i,j)*uuu(i,j)
//         rov(i,j)=rho(i,j)*vvv(i,j)
//         roe(i,j)=rho(i,j)*eee(i,j)
//         scp(i,j)=1.
//         x=x+dlx
//      enddo
//      y=y+dly
//   enddo
