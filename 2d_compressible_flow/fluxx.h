//

void fluxx1(ACC<double> fro, const ACC<double> tb1, const ACC<double> tb2) {
    fro(0,0)=-tb1(0,0)-tb2(0,0);
}

void fluxx2(ACC<double> tb1, const ACC<double> rou, const ACC<double> uuu) {
    tb1(0,0)=rou(0,0)*uuu(0,0);
}

void fluxx3(ACC<double> &tba, const ACC<double> &tb6, const ACC<double> &tb7, const ACC<double> &tb9) {
    tba(0,0)=xmu*(qtt*tb6(0,0)+tb7(0,0)+utt*tb9(0,0));
}

void fluxx4(ACC<double> &fru, const ACC<double> &tb3, const ACC<double> &tb4, const ACC<double> &tb5, const ACC<double> &tba, const ACC<double> &eps, const ACC<double> &uuu) {
    fru(0,0)=-tb3(0,0)-tb4(0,0)-tb5(0,0)+tba(0,0)-(eps(0,0)/eta)*uuu(0,0);
}