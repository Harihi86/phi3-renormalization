// Exact transport around the full angular moment's t=v ridge.
#include "g371_fullangle_functions_v1.hpp"
Real ridge_host_value(const DoubleHostGraph&g,bool positive,bool tail,const std::array<Real,6>&p){
    auto pp=p;Real jac=1;
    if(tail){Real den=sq(p[0])+sq(1-p[0]);pp[0]=sq(p[0])/den;jac=2*p[0]*(1-p[0])/sq(den);}
    Real total=sq(pp[0]/(1-pp[0])),root=rt(1+2*total),scale=root/total,Y;
    arb_asinh(Y.v,(total/root).v,precision);Real eta=(2*p[1]-1)*Y,sh,ch;
    arb_sinh_cosh(sh.v,ch.v,eta.v,precision);pp[1]=(1+scale*sh)/2;
    pos(pp[1],"positive transported first fraction");pos(1-pp[1],"positive transported second fraction");
    return jac*scale*Y*ch*fullangle_host_value(g,positive,true,true,pp);
}
