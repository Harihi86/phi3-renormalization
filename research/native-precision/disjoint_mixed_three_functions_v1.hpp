// Exact original mixed child product after three complete loop integrations.
struct MixedThreeGraph{std::string id,sha;std::array<std::array<int,3>,3>host,triangle;std::array<std::array<int,3>,2>bubbles;};
#include "disjoint_mixed_three_data.hpp"
struct OneHostJet{Real value,constant;std::array<Real,2>linear;std::array<std::array<Real,2>,2>square;};
OneHostJet multiply_one_host(const OneHostJet&a,const OneHostJet&b){
    OneHostJet c;c.value=a.value*b.value;c.constant=a.constant*b.value+a.value*b.constant;
    for(int i=0;i<2;++i){c.linear[i]=a.linear[i]*b.value+a.value*b.linear[i];for(int j=0;j<2;++j)c.square[i][j]=a.value*b.square[i][j]+b.value*a.square[i][j]+(a.linear[i]*b.linear[j]+b.linear[i]*a.linear[j])/2;}return c;
}
std::array<Real,3> mixed_bubble_jets(const Real&t){
    pos(t,"positive complete bubble argument");Real E=edge_log_full(t),a;arb_asinh(a.v,(rt(t)/2).v,precision);
    Real j0=4*a/rt(t*(t+4)),B=(t+4)*E/6-t/9,B1=E/6+(t+4)*(1-j0)/(6*t)-Real(1)/9,B2=Real(1)/(6*t)-(1-j0)/sq(t),D=1+t;
    return {B/D,B1/D-B/sq(D),B2/D-2*B1/sq(D)+2*B/power3(D)};
}
Real mixed_bubble_value(const Real&t){pos(t,"positive shifted bubble invariant");return ((t+4)*edge_log_full(t)/6-t/9)/(1+t);}
Real mixed_three_value(const MixedThreeGraph&g,bool positive,bool control,const std::array<Real,3>&p){
    Real s=p[0],z=p[1],theta=p[2],t=sq(s/(1-s)),A=1+t,co,si;
    arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real rad=rt(t*z),sqrt3=rt(Real(3)),T0=-edge_log_full(t)/2,total;
    auto bj=mixed_bubble_jets(t);Real J=bj[0],P0=sq(J)*T0;
    for(long sign:{-1L,1L}){
        Real qx=Real(sign)*rad*co,qy=Real(sign)*rad*si,q1=qx,q2=(-qx+sqrt3*qy)/2;
        auto norm=[&](const std::array<int,3>&r){
            if(positive)return sq(Real(r[0])*qx+Real(r[1])-Real(r[2])/2)+sq(Real(r[0])*qy+sqrt3*Real(r[2])/2)+Real(r[0]*r[0])*t*(1-z);
            return Real(r[0]*r[0])*t+2*Real(r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real lr;for(const auto&r:g.host){Real prop=1+norm(r);pos(prop,"positive mixed host denominator");lr=lr-lp((prop-A)/A);}
        Real d0=mixed_bubble_value(norm(g.bubbles[0]))-J,d1=mixed_bubble_value(norm(g.bubbles[1]))-J;
        Real dt=-triangle_log_full(norm(g.triangle[0]),norm(g.triangle[1]))-T0;
        Real difference=T0*(J*(d0+d1)+d0*d1)+dt*(J+d0)*(J+d1);
        total=total+(em1(lr)*P0+real_exp(lr)*difference)/2;
    }
    Real common=2/(s*(1-s))*power3(t/A),answer=common*(1-z)*total;
    if(control){
        auto extnorm=[](const std::array<int,3>&r){return Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        auto lv=[](const std::array<int,3>&r){return std::array<Real,2>{Real(2*r[0]*r[1]),Real(2*r[0]*r[2])};};
        OneHostJet product;product.value=1;
        for(const auto&r:g.host){OneHostJet one;one.value=1;auto l=lv(r);one.constant=-extnorm(r)/A;
            for(int i=0;i<2;++i){one.linear[i]=-l[i]/A;for(int j=0;j<2;++j)one.square[i][j]=l[i]*l[j]/sq(A);}product=multiply_one_host(product,one);}
        for(const auto&r:g.bubbles){OneHostJet one;one.value=J;auto l=lv(r);one.constant=bj[1]*extnorm(r);
            for(int i=0;i<2;++i){one.linear[i]=bj[1]*l[i];for(int j=0;j<2;++j)one.square[i][j]=bj[2]*l[i]*l[j]/2;}product=multiply_one_host(product,one);}
        Real as;arb_asinh(as.v,(rt(t)/2).v,precision);Real j0=4*as/rt(t*(t+4)),j2=2*(1+j0)/(t+4),ell=-(1-j0)/(2*t),Cp=ell/2,Cf=-(j0-3*(1-j0)/t)/12,I2=(1-2*j0+j2)/sq(t);
        auto a=lv(g.triangle[0]),b=lv(g.triangle[1]);for(int i=0;i<2;++i){a[i]=a[i]/2;b[i]=b[i]/2;}
        OneHostJet triangle;triangle.value=T0;triangle.constant=Cp*(extnorm(g.triangle[0])+extnorm(g.triangle[1]))+Cf*extnorm(g.triangle[2]);
        for(int i=0;i<2;++i){triangle.linear[i]=ell*(a[i]+b[i]);for(int j=0;j<2;++j)triangle.square[i][j]=I2*(a[i]*a[j]+b[i]*b[j]+(a[i]*b[j]+b[i]*a[j])/2)/3;}
        product=multiply_one_host(product,triangle);std::array<Real,2>q={rad*co,rad*(-co+sqrt3*si)/2};Real point=product.constant;
        for(int i=0;i<2;++i)for(int j=0;j<2;++j)point=point+product.square[i][j]*q[i]*q[j];
        Real mean=product.constant+t*(product.square[0][0]+product.square[1][1]-(product.square[0][1]+product.square[1][0])/2)/6;
        answer=answer+common*(mean/2-(1-z)*point);
    }
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite mixed three-cube forest required");return answer;
}
MixedThreeGraph mixed_three_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original recipe required");auto g=load_mixed_three(id);if(sha!=g.sha)throw std::invalid_argument("original mixed graph identity required");return g;}
