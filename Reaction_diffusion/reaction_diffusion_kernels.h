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


void initcond_stencil(ACC<double> &A, const int *idx) {
    double a = 0;
    for (int m=0; m < 5; m++) {
    a = a + 10*pow(-1, m)*4/(2*m+1)/(2*m+1)/pi/pi*sin((2*m+1)*pi*(idx[0]*dx));
    
    }
    A(0, 0) = a;
    
}

void laplacian(ACC<double> &A) {
    
}