#include "kernel.hpp"
#include <iostream>
using namespace nested_e2;
void show(const char*name,const J2&a){std::cout<<name;for(auto&c:a.c){std::cout<<" ";arb_printn(c.v,12,0);}std::cout<<'\n';}
J2 rem(const J2&t,const J2&p){
 J2 D=t+1,difference=D-p;show("D",D);show("difference",difference);
 J2 sq=square(difference);show("square",sq);J2 arg=sq+4*p;show("arg",arg);
 J2 radical=root(arg);show("radical",radical);J2 delta;
 if(arb_is_positive(difference.c[0].v)){show("rational denominator",difference+radical);delta=2*p/(difference+radical);}else delta=(p-D+radical)/2;
 show("delta",delta);J2 E=(D+p+radical)/2;show("E",E);show("D3",cube(D));show("E3",cube(E));
 J2 r=p*delta*(3*square(D)+6*D*delta+(t+3)*square(delta))/3/cube(D)/cube(E);show("remainder",r);return r;
}
int main(){precision=192;Geometry g{234,8,225,8};Real xa=dyadic(g.ix,g.lx,24),xb=dyadic(g.ix+1,g.lx,24),ya=dyadic(g.iy,g.ly,25),yb=dyadic(g.iy+1,g.ly,25),xx,yy;arb_union(xx.v,xa.v,xb.v,precision);arb_union(yy.v,ya.v,yb.v,precision);
 for(int axis=0;axis<2;++axis){std::cout<<"AXIS "<<axis<<'\n';J2 x(xx),y(yy);if(axis==0)x.c[1]=1;else y.c[1]=1;
 J2 t=expj(J2(arb_cw::log(Real(2)))*x),u=expj(J2(arb_cw::log(Real(2)))*y);show("t",t);show("u",u);show("bubble",bubble(u));std::cout<<"OUTER\n";rem(t,1);std::cout<<"INNER\n";rem(u,t);show("density",density(x,y));}
}
