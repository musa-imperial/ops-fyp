//kernels

void set_zero(ACC<double> &A) {
    A(0, 0) = 0.0;
}

void copy(ACC<double> &A, const ACC<double> &Anew) {
    A(0, 0) = Anew(0, 0);
}

void initial_condition(ACC<double> &A, const int *idx) {
     
    A(0,0) = 5*sin(M_PI*((idx[0]+1)*dx)/Lx)*sin(M_PI*((idx[1]+1)*dy)/Ly);
}

//top, bottom, left and right boundary condition kernels
void boundary_condition(ACC<double> &Anew) {
    Anew(0,0) = 0;
}

void interior_stencil_u(const ACC<double> &A, ACC<double> &Anew) {
    //Anew(0, 0) = A(0, 0)+hmudt*(A(1, 0)+A(-1, 0)+A(0, 1)+A(0, -1)-4*A(0, 0));
    Anew(0, 0) = a0*A(0, 0)+hmudt*(A(1, 0)+A(-1, 0)+A(0, 1)+A(0, -1));
}

void solution_check(const ACC<double> &A, const int *idx, double *l2error, double *max_error) {
    double u = 5*exp(-mu*M_PI*M_PI*(1/Lx/Lx+1/Ly/Ly)*T)*sin(M_PI/Lx*(idx[0]+1)*dx)*sin(M_PI/Ly*(idx[1]+1)*dy);
    *l2error = *l2error + sqrt(A(0,0)*A(0,0)-u*u);
    *max_error = fmax(*max_error, fabs((A(0,0)-u)/u));
    //std::cout << (idx[0]+1)*dx << ", "<< (idx[1]+1)*dy << ", " << (A(0,0)) << ", " << u << ", " << fabs((A(0,0)-u))/u << std::endl;
}
