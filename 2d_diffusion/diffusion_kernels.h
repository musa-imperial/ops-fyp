//kernels

void set_zero(ACC<double> &A) {
    A(0, 0) = 0.0;
}

void copy(ACC<double> &A, const ACC<double> &Anew) {
    A(0, 0) = Anew(0, 0);
}

void u_initcond_stencil(ACC<double> &A, const int *idx) {

    if ((Lx/4 < idx[0] && Lx/4+a < idx[0]) && ((Lx/4 < idx[0] && Lx/4+a < idx[0])) && (Ly/4 < idx[1] && Ly/4+a < idx[1]) && ((Ly/4 < idx[1] && Ly/4+a < idx[1]))) {
        A(0, 0) = 1;
    }
    else {
        A(0, 0) = 0;
    }
    
}


void bottomleft_u(const ACC<double> &A, ACC<double> &Anew) {
    
     Anew(0, 0) = A(0, 0)+(hnudt*(A(1, 0) + A(0, 1) - 2 * A(0, 0)));
    
}
void topleft_u(const ACC<double> &A, ACC<double> &Anew) {
    
     Anew(0, 0) = A(0, 0)+(hnudt*(A(1, 0) + A(0, -1) - 2 * A(0, 0)));

}
void bottomright_u(const ACC<double> &A, ACC<double> &Anew) {
    
     Anew(0, 0) = A(0, 0)+(hnudt*(A(-1, 0) + A(0, 1) - 2 * A(0, 0)));

}

void topright_u(const ACC<double> &A, ACC<double> &Anew) {

     Anew(0, 0) = A(0, 0)+(hnudt*(A(-1, 0) + A(0, -1) - 2 * A(0, 0)));
}


//top, bottom, left and right boundary condition kernels
void left_bndcon_u(const ACC<double> &A, ACC<double> &Anew) {
    Anew(0, 0) = A(0, 0)+hnudt*(A(1, 0) + A(0, 1) + A(0, -1) - 3 * A(0, 0));
}

void right_bndcon_u(const ACC<double> &A, ACC<double> &Anew) {
    Anew(0, 0) = A(0, 0)+hnudt*(A(-1, 0) + A(0, 1) + A(0, -1) - 3 * A(0, 0));
}

////////
void top_bndcon_u(const ACC<double> &A, ACC<double> &Anew) {
    Anew(0, 0) = A(0, 0)+hnudt*(A(1, 0) + A(-1, 0) + A(0, -1) - 3 * A(0, 0));
}


void bottom_bndcon_u(const ACC<double> &A, ACC<double> &Anew) {
    Anew(0, 0) = A(0, 0)+hnudt*(A(1, 0) + A(-1, 0) + A(0, 1) - 3 * A(0, 0));
}


void interior_stencil_u(const ACC<double> &A, ACC<double> &Anew) {
    Anew(0, 0) = A(0, 0)+hnudt*(A(1, 0)+A(-1, 0)+A(0, 1)+A(0, -1)-4*A(0, 0));
}

