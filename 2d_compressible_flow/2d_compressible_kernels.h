void set_zero(ACC<double> &A) {
    A(0, 0) = 0.0;
}

void copy(ACC<double> &A, const ACC<double> &Anew) {
    A(0, 0) = Anew(0, 0);
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
    if ((idx[0]-xlx/2)*(idx[0]-xlx/2) + (idx[1]-yly/2)*(idx[1]-yly/2) < radius*radius) {
        eps(0, 0) = 1;
    }
    else {
        eps(0, 0) = 0;
    }
}

