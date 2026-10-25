// Order four is a zero-integral angular control, never a new local subtraction.
#include "nested_two_bubble_functions_v1.hpp"
using J4=std::array<Real,5>;
J4 j4constant(const Real&a){J4 x;x[0]=a;return x;}
J4 j4add(const J4&a,const J4&b){J4 c;for(int i=0;i<5;++i)c[i]=a[i]+b[i];return c;}
J4 j4scale(const J4&a,const Real&b){J4 c;for(int i=0;i<5;++i)c[i]=a[i]*b;return c;}
J4 j4mul(const J4&a,const J4&b){J4 c;for(int i=0;i<5;++i)for(int j=0;j<=i;++j)c[i]=c[i]+a[j]*b[i-j];return c;}
J4 j4inv(const J4&a){J4 c;c[0]=1/a[0];for(int i=1;i<5;++i){Real s;for(int j=1;j<=i;++j)s=s+a[j]*c[i-j];c[i]=-s/a[0];}return c;}
J4 j4sqrt(const J4&a){J4 c;c[0]=rt(a[0]);for(int i=1;i<5;++i){Real s;for(int j=1;j<i;++j)s=s+c[j]*c[i-j];c[i]=(a[i]-s)/(2*c[0]);}return c;}
J4 j4asinh(const J4&a){J4 c,der;arb_asinh(c[0].v,a[0].v,precision);for(int i=0;i<4;++i)der[i]=Real(i+1)*a[i+1];
    J4 f=j4mul(der,j4inv(j4sqrt(j4add(j4constant(Real(1)),j4mul(a,a)))));for(int i=1;i<5;++i)c[i]=f[i-1]/Real(i);return c;}
J4 j4bubble(const J4&a){J4 plus=j4add(a,j4constant(Real(4))),root=j4sqrt(a);
    return j4add(j4add(j4scale(j4mul(j4mul(plus,j4sqrt(plus)),j4mul(j4inv(root),j4asinh(j4scale(root,Real(1)/2)))),Real(1)/3),j4scale(a,Real(-4)/9)),j4constant(Real(-4)/3));}
J4 j4compose(const J4&a,const J4&d){J4 c=j4constant(a[4]);for(int i=3;i>=0;--i)c=j4add(j4mul(c,d),j4constant(a[i]));return c;}
Real two_bubble_quartic_value(const TwoBubbleGraph&g,bool positive,const std::array<Real,5>&p){
    Real result=two_bubble_value(g,positive,true,p),s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],t=sq(s/(1-s)),A=1+t;
    Real den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den),co,si,sqrt3=rt(Real(3));arb_sin_cos_pi(si.v,co.v,theta.v,precision);
    std::array<std::array<int,3>,3>parent={g.parent[0],g.parent[1],g.parent[0]};for(int j=0;j<3;++j)parent[2][j]-=g.parent[1][j];
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};Real V0,Cp;std::array<Real,2>pv;
    for(int j=0;j<3;++j){const auto&r=parent[j];V0=V0+w[j]*Real(r[0]*r[0])*t;Cp=Cp+w[j]*Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);for(int k=0;k<2;++k)pv[k]=pv[k]+w[j]*Real(2*r[0]*r[k+1]);}
    J4 argument=j4constant((1-x+V0)/x);argument[1]=1/x;J4 denominator=j4constant(1+V0);denominator[1]=1;
    J4 deriv=j4scale(j4mul(j4bubble(argument),j4mul(j4inv(denominator),j4inv(denominator))),-x),K;
    K[0]=-(bubble_A((1-x+V0)/x)-bubble_A((1-x)/x));for(int n=1;n<5;++n)K[n]=deriv[n-1]/Real(n);
    J4 tj=j4constant(t);tj[1]=1;J4 J=j4mul(j4bubble(tj),j4inv(j4add(tj,j4constant(Real(1)))));
    auto coefficient=[&](const Real&qx,const Real&qy){
        std::array<Real,2>q={qx,(-qx+sqrt3*qy)/2};J4 product=j4constant(Real(1));
        for(const auto&r:g.host){J4 d=j4constant(Real(1));d[1]=2*Real(r[0])*(Real(r[1])*q[0]+Real(r[2])*q[1])/A;d[2]=Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2])/A;product=j4mul(product,j4inv(d));}
        J4 pd;pd[1]=pv[0]*q[0]+pv[1]*q[1];pd[2]=Cp;product=j4mul(product,j4compose(K,pd));
        const auto&r=g.bubble;J4 bd;bd[1]=2*Real(r[0])*(Real(r[1])*q[0]+Real(r[2])*q[1]);bd[2]=Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);product=j4mul(product,j4compose(J,bd));return product[4];};
    Real point=coefficient(rt(t*z)*co,rt(t*z)*si),mean=coefficient(Real(0),Real(0))/3;
    for(int i=0;i<8;++i){Real c,h;arb_sin_cos_pi(h.v,c.v,(Real(i)/4).v,precision);mean=mean+coefficient(rt(t/2)*c,rt(t/2)*h)/12;}
    Real common=2/(s*(1-s))*power3(t/A)*(1-x)*jac;
    result=result+common*(mean/2-(1-z)*point);if(!arb_is_finite(result.v))throw std::runtime_error("finite quartic angular control required");return result;
}
Real two_bubble_quartic_mapped(const TwoBubbleGraph&g,bool positive,bool matched,const std::array<Real,5>&p){
    if(!matched)return two_bubble_quartic_value(g,positive,p);
    Real x=sq(p[3])/(sq(p[3])+sq(1-p[3])),y=p[4],s=p[0];std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    std::array<int,3>loops={g.parent[0][0],g.parent[1][0],g.parent[0][0]-g.parent[1][0]};Real ell;for(int j=0;j<3;++j)ell=ell+w[j]*Real(loops[j]*loops[j]);
    pos(ell,"positive matched parent variance");Real h=1/rt(ell),den=1-s+h*s;auto pp=p;pp[0]=h*s/den;return h/sq(den)*two_bubble_quartic_value(g,positive,pp);
}
