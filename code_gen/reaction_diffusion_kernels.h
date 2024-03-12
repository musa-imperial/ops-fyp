//kernels

void set_zero(ACC<double> &A) {
    A(0, 0) = 0.0;
}

void copy(ACC<double> &A, const ACC<double> &Anew) {
    A(0, 0) = Anew(0, 0);
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


void bottom_left_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu1dt*(A(1, 0) + A(0, 1) - 2 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt);

}

void bottom_left_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu2dt*(A(1, 0) + A(0, 1) - 2 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt);

}

void top_left_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu1dt*(A(1, 0) + A(0, -1) - 2 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt);

}

void top_left_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu2dt*(A(1, 0) + A(0, -1) - 2 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt);

}

void bottom_right_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu1dt*(A(-1, 0) + A(0, 1) - 2 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt);

}

void bottom_right_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu2dt*(A(-1, 0) + A(0, 1) - 2 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt);

}

void top_right_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu1dt*(A(-1, 0) + A(0, -1) - 2 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt);

}

void top_right_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    
     Anew(0, 0) = A(0, 0)+(hmu2dt*(A(-1, 0) + A(0, -1) - 2 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt);

}

//top, bottom, left and right boundary condition kernels
void left_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu1dt*(A(1, 0) + A(0, 1) + A(0, -1) - 3 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt;
}

void left_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu2dt*(A(1, 0) + A(0, 1) + A(0, -1) - 3 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt;
}

void right_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu1dt*(A(-1, 0) + A(0, 1) + A(0, -1) - 3 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt;
}

void right_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu2dt*(A(-1, 0) + A(0, 1) + A(0, -1) - 3 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt;
}
////////
void top_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu1dt*(A(1, 0) + A(-1, 0) + A(0, -1) - 3 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt;
}

void top_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu2dt*(A(1, 0) + A(-1, 0) + A(0, -1) - 3 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt;
}

void bottom_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu1dt*(A(1, 0) + A(-1, 0) + A(0, 1) - 3 * A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt;
}

void bottom_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu2dt*(A(1, 0) + A(-1, 0) + A(0, 1) - 3 * A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt;
}

void interior_stencil_u(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu1dt*(A(1, 0)+A(-1, 0)+A(0, 1)+A(0, -1)-4*A(0, 0))+(eps * A(0, 0) * (1 - A(0, 0)) * (A(0, 0) - (B(0, 0) + b) * div_a)) * dt;
}

void interior_stencil_v(const ACC<double> &A, ACC<double> &Anew, const ACC<double> &B) {
    Anew(0, 0) = A(0, 0)+hmu2dt*(A(1, 0)+A(-1, 0)+A(0, 1)+A(0, -1)-4*A(0, 0))+(B(0, 0) * B(0, 0) * B(0, 0) - A(0, 0)) * dt;
}