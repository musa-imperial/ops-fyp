

void derix(const ACC<double> &phi, ACC<double> &dfi, const double *udx ) {
    //{0,0)
    //{1,0, -1,0}
    dfi(0,0)=(*udx)*(phi(1,0)-phi(-1,0));
}

void deriy(const ACC<double> &phi, ACC<double> &dfi, const double *udy ) {
    //{0,0)
    //{0,1, 0,-1}
    dfi(0,0)=(*udy)*(phi(0,1)-phi(0,-1));
}

// void derix_periodic(const ACC<double> &phi, ACC<double> &dfi ) {
//     //{0,0)
//     //{1,0, -1,0}
//     dfi(0,0)=udx*(phi(1,0)-phi(-1,0));
// }


// void deriy(const ACC<double> &phi, ACC<double> &dfi ) {
//     //{0,0)
//     //{0,1, 0,-1}
//     dfi(0,0)=udy*(phi(0,1)-phi(0,-1));
// }

// void deriy_periodic(const ACC<double> &phi, ACC<double> &dfi ) {
//     //{0,0)
//     //{0,1, 0,-1}
//     dfi(0,0)=udy*(phi(0,1)-phi(0,-1));
// }

void derxx(const ACC<double> &phi, ACC<double> &dfi, const double *udx ) {
    //{0,0)
    //{0,0, 1,0, -1,0}
    dfi(0,0)=(*udx)*(phi(1,0)-(phi(0,0)+phi(0,0))+phi(-1,0));
}

void deryy(const ACC<double> &phi, ACC<double> &dfi, const double *udx ) {
    //{0,0)
    //{0,0, 1,0, -1,0}
    dfi(0,0)=(*udx)*(phi(0,1)-(phi(0,0)+phi(0,0))+phi(0,-1));
}