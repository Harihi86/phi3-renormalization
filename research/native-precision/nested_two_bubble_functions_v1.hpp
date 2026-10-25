// All forty original forest terms after exact elimination of three child loops.
struct TwoBubbleGraph{std::string id,sha;std::array<std::array<int,3>,2>parent;std::array<std::array<int,3>,3>host;std::array<int,3>bubble;};
#include "nested_two_bubble_data.hpp"
Real two_bubble_value(const TwoBubbleGraph&g,bool positive,bool control,const std::array<Real,5>&p){
    Real s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],t=sq(s/(1-s)),A=1+t;
    Real den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den),co,si;
    arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real radial=rt(t*z),sqrt3=rt(Real(3));
    std::array<std::array<int,3>,3>parent={g.parent[0],g.parent[1],g.parent[0]};for(int j=0;j<3;++j)parent[2][j]-=g.parent[1][j];
    if(g.bubble[0]!=1)throw std::invalid_argument("normalized disjoint bubble route required");
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    Real V0;for(int j=0;j<3;++j)V0=V0+w[j]*Real(parent[j][0]*parent[j][0])*t;
    Real a0=bubble_A((1-x+V0)/x),K0=-(a0-bubble_A((1-x)/x)),bt=bubble_value(t),J0=bt/A,total;
    for(long sign:{-1L,1L}){
        Real qx=Real(sign)*radial*co,qy=Real(sign)*radial*si,q1=qx,q2=(-qx+sqrt3*qy)/2;
        auto norm=[&](const std::array<int,3>&r){
            if(positive)return sq(Real(r[0])*qx+Real(r[1])-Real(r[2])/2)+sq(Real(r[0])*qy+sqrt3*Real(r[2])/2)+Real(r[0]*r[0])*t*(1-z);
            return Real(r[0]*r[0])*t+2*Real(r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real V,lr;for(int j=0;j<3;++j){V=V+w[j]*norm(parent[j]);Real host=1+norm(g.host[j]);pos(host,"positive distinct host propagator");lr=lr-lp((host-A)/A);}
        Real bnorm=norm(g.bubble);pos(bnorm,"positive massive bubble invariant in open cube");
        Real DK=-(bubble_A((1-x+V)/x)-a0),DJ=bubble_value(bnorm)/(1+bnorm)-J0,ratio=real_exp(lr);
        total=total+(em1(lr)*K0*J0+ratio*(DK*J0+K0*DJ+DK*DJ))/2;
    }
    Real common=2/(s*(1-s))*power3(t/A)*(1-x)*jac,answer=common*(1-z)*total;
    if(control){
        std::array<Real,2>q={radial*co,radial*(-co+sqrt3*si)/2},pv,hv;
        Real Cp,hsquare,hconstant,hpoint;
        for(int j=0;j<3;++j){const auto&a=parent[j];const auto&h=g.host[j];
            Cp=Cp+w[j]*Real(a[1]*a[1]+a[2]*a[2]-a[1]*a[2]);
            std::array<Real,2>one={Real(2*h[0]*h[1]),Real(2*h[0]*h[2])};hsquare=hsquare+norm_gram(one);
            hconstant=hconstant+Real(h[1]*h[1]+h[2]*h[2]-h[1]*h[2]);hpoint=hpoint+sq(one[0]*q[0]+one[1]*q[1]);
            for(int k=0;k<2;++k){pv[k]=pv[k]+w[j]*Real(2*a[0]*a[k+1]);hv[k]=hv[k]+one[k];}}
        std::array<Real,2>rv={-hv[0]/A,-hv[1]/A},bv={Real(2*g.bubble[1]),Real(2*g.bubble[2])};
        Real Cb=Real(g.bubble[1]*g.bubble[1]+g.bubble[2]*g.bubble[2]-g.bubble[1]*g.bubble[2]);
        Real Lp=q[0]*pv[0]+q[1]*pv[1],Lb=q[0]*bv[0]+q[1]*bv[1],R1=q[0]*rv[0]+q[1]*rv[1];
        Real R2=(sq(q[0]*hv[0]+q[1]*hv[1])+hpoint)/(2*sq(A))-hconstant/A;
        Real meanR2=t*(norm_gram(hv)+hsquare)/(12*sq(A))-hconstant/A;
        Real arg=(1-x+V0)/x,B=bubble_value(arg),Kp=-x*B/sq(1+V0),Kpp=2*x*B/power3(1+V0)-bubble_derivative(arg)/sq(1+V0);
        Real hp;arb_asinh(hp.v,(rt(t)/2).v,precision);Real edge=4*hp/rt(t*(t+4));
        Real bp=bubble_derivative(t),bpp=1/(6*t)-(1-edge)/sq(t);
        Real Jp=bp/A-bt/sq(A),Jpp=bpp/A-2*bp/sq(A)+2*bt/power3(A);
        Real point=R2*K0*J0+R1*(Kp*Lp*J0+K0*Jp*Lb)+(Kp*Cp+Kpp*sq(Lp)/2)*J0+Kp*Jp*Lp*Lb+K0*(Jp*Cb+Jpp*sq(Lb)/2);
        Real mean=meanR2*K0*J0+t/6*(dot_gram(rv,pv)*Kp*J0+dot_gram(rv,bv)*K0*Jp)+(Kp*Cp+Kpp*t*norm_gram(pv)/12)*J0+Kp*Jp*t*dot_gram(pv,bv)/6+K0*(Jp*Cb+Jpp*t*norm_gram(bv)/12);
        answer=answer+common*(mean/2-(1-z)*point);
    }
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite complete two-bubble five-cube required");return answer;
}
Real two_bubble_mapped(const TwoBubbleGraph&g,bool positive,bool matched,bool control,const std::array<Real,5>&p){
    if(!matched)return two_bubble_value(g,positive,control,p);
    Real x=sq(p[3])/(sq(p[3])+sq(1-p[3])),y=p[4],s=p[0];
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    std::array<int,3>loops={g.parent[0][0],g.parent[1][0],g.parent[0][0]-g.parent[1][0]};Real ell;
    for(int j=0;j<3;++j)ell=ell+w[j]*Real(loops[j]*loops[j]);pos(ell,"positive exact parent variance");
    Real h=1/rt(ell),den=1-s+h*s;auto pp=p;pp[0]=h*s/den;
    return h/sq(den)*two_bubble_value(g,positive,control,pp);
}
TwoBubbleGraph two_bubble_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original two-bubble graph recipe");auto g=load_two_bubble(id);if(g.sha!=sha)throw std::invalid_argument("original two-bubble graph hash required");return g;}
