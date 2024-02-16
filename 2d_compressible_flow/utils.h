void param(double *xlx, double *yly, double *xmu, double *xba, double *gma,
                  double *chp, double *roi, double *cci, double *d, double *tpi,
                  double *chv, double *uu0) {
    double ren, pdl, mach;

    ren=200. //REYNOLDS NUMBER
    mach=0.2 //MACH NUMBER
    pdl=0.7  //PRANDTL NUMBER
    roi=1.
    cci=1.   //SOUND SPEED
    d=1.
    chp=1.
    gma=1.4
        
    chv=chp/gma
    xlx=4.*d     //DOMAIN SIZE X DIRECTION
    yly=4.*d     //DOMAIN SIZE Y DIRECTION
    uu0=mach*cci
    xmu=roi*uu0*d/ren
    xba=xmu*chp/pdl
    tpi=cci**2/(chp*(gma-1))

}