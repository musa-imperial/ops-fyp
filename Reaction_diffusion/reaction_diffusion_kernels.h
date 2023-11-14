#define laplacianu(a, b, c, d) ((a+b+c-3*d))
#define laplacianv(a, b, c, d) ((a+b+c-3*d))
#define laplacian_corner(a, b, c) ((a+b-2*c))
#define laplacianvc(a, b, c) ((a+b-2*c))
#define f1(u, v) (epsilon*a*(1-u)*(u-(v+b)*div_a))*dt
#define f2(v, u) (u*u*u-v)*dt

void set_zero(ACC<double> &A) {
    A(0, 0) = 0.0;
}

void copy(ACC<double> &A, const ACC<double> &Anew) {
    A(0, 0) = Anew(0, 0);
}

void left_bndcon(ACC<double> &A, const int *idx) {
    A(0, 0) = 0;
}

void right_bndcon(ACC<double> &A, const int *idx) {
    A(0, 0) = A(-1, 0);
}

void apply_stencil(const ACC<double> &A, ACC<double> &Anew) {
    Anew(0, 0) = C1*A(1, 0)+C2*A(0, 0)+C3*A(-1, 0);
}


void u_initcond_stencil(ACC<double> &A, const int *idx) {

    if (idx[1] > Ly/2) {
        A(0, 0) = 1;
    }
    else {
        A(0, 0) = 0;
    }
    
}

void v_initcond_stencil(ACC<double> &A, const int *idx) {

    if (idx[0] < Lx/2) {
        A(0, 0) = a/2;
    }
    else {
        A(0, 0) = 0;
    }
    
}


void bottomleft_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu1dt*laplacian_corner(A(1, 0), A(0, 1), A(0, 0))+f1(A(0, 0), B(0, 0)))+ A(0, 0);

}

void bottomleft_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu2dt*laplacian_corner(A(1, 0), A(0, 1), A(0, 0))+f2(A(0, 0), B(0, 0)))+ A(0, 0);

}

void topleft_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu1dt*laplacian_corner(A(1, 0), A(0, -1), A(0, 0))+f1(A(0, 0), B(0, 0)))+ A(0, 0);

}

void topleft_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu2dt*laplacian_corner(A(1, 0), A(0, -1), A(0, 0))+f2(A(0, 0), B(0, 0)))+ A(0, 0);

}

void bottomright_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu1dt*laplacian_corner(A(-1, 0), A(0, 1), A(0, 0))+f1(A(0, 0), B(0, 0)))+ A(0, 0);

}

void bottomright_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu2dt*laplacian_corner(A(-1, 0), A(0, 1), A(0, 0))+f2(A(0, 0), B(0, 0)))+ A(0, 0);

}

void topright_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu1dt*laplacian_corner(A(-1, 0), A(1, 0), A(0, 0))+f1(A(0, 0), B(0, 0)))+ A(0, 0);

}

void topright_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = (hmu2dt*laplacian_corner(A(-1, 0), A(1, 0), A(0, 0))+f2(A(0, 0), B(0, 0)))+ A(0, 0);

}