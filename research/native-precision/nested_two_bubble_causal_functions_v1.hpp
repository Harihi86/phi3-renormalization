// Parameterize the actual parent variance before applying the matched radial map.
#include "nested_two_bubble_quartic_functions_v1.hpp"
Real two_bubble_causal_value(const TwoBubbleGraph&g,bool positive,bool balanced,const std::array<Real,5>&p){
    Real u=p[3],v=p[4],du=sq(u)+sq(1-u),xi=sq(u)/du,jxi=2*u*(1-u)/sq(du),w=v,jw=1;
    if(balanced){Real dv=sq(v)+sq(1-v);w=sq(v)/dv;jw=2*v*(1-v)/sq(dv);}
    Real x=(1-xi)*w,y=xi/(1-x),sx=rt(x),cx=rt(1-x),old=sx/(sx+cx);
    Real oldden=sq(old)+sq(1-old),oldjac=2*old*(1-old)/sq(oldden);
    auto pp=p;pp[3]=old;pp[4]=y;pos(y,"positive causal parent fraction");pos(1-y,"strict causal parent fraction");
    Real jac=jxi*jw*(1-xi)/((1-x)*oldjac);
    return jac*two_bubble_quartic_mapped(g,positive,true,pp);
}
