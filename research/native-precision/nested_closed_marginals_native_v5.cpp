// Complete original forest after exact parent and triangle scalar marginals.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <cmath>
#include <algorithm>
namespace nested_closed {
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::Z;
using fourier_e2_root::readq;using fourier_e2_root::interval;
struct Graph{std::string id,sha;std::array<std::array<int,3>,2> parent;std::array<std::array<int,3>,3> triangle,host;};
#include "nested_mixed_data.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real lp(const Real&a){Real r;arb_log1p(r.v,a.v,precision);return r;}
Real real_log(const Real&a){Real r;arb_log(r.v,a.v,precision);return r;}
Real real_exp(const Real&a){Real r;arb_exp(r.v,a.v,precision);return r;}
Real em1(const Real&a){Real r;arb_expm1(r.v,a.v,precision);return r;}
Real power3(const Real&a){return a*sq(a);}
void pos(const Real&a,const char*msg){if(!arb_is_positive(a.v))throw std::runtime_error(msg);}
Z clog(const Z&a){Z r;acb_log(r.v,a.v,precision);return r;}
Z dilog(const Z&a){Z r;acb_polylog_si(r.v,2,a.v,precision);return r;}
Real atan_log(const Real&q){
    Z a;arb_sqrt_ui(acb_imagref(a.v),3,precision);Z sum;
    for(long b:{-1L,1L}){Z zz=-a/(Z(b)-a);Real u=1+Real(b)*q;pos(u,"positive marginal logarithm argument");
        sum=sum+(clog(Z(1)-a/Z(b))*Z(real_log(u))-dilog(zz*Z(u))+dilog(zz))/Z(b);}
    Real r;arb_set(r.v,acb_imagref(sum.v));return r/2;
}
Real bubble_A(const Real&t){
    if(arb_is_zero(t.v))return Real(0);pos(t,"nonnegative bubble primitive argument");
    Real q=rt(t/(t+4)),h,angle;pos(1-q,"bubble primitive finite q");
    arb_asinh(h.v,(rt(t)/2).v,precision);arb_atan(angle.v,(rt(Real(3))*q).v,precision);
    Real value=sq(h)/3+2*rt(Real(3))*(h*angle-atan_log(q))+4*q*h/(1+3*sq(q))-Real(17)/18*lp(t)-Real(32)/9*sq(q)/(1+3*sq(q));
    if(!arb_is_finite(value.v))throw std::runtime_error("finite bubble primitive required");return value;
}
Real linear_mean(const Real&z){
    if(arb_is_zero(z.v))return Real(0);pos(1+z,"positive linear logarithm");
    Real abs;arb_abs(abs.v,z.v);
    if(arb_lt(abs.v,(Real(1)/2).v)){
        // P(z)=sum_{j>=1}(-1)^(j+1) z^j/[j(j+1)].
        // For r=|z|<1/2 the omitted tail is <=2*r^(N+1).
        mag_t radius,tail;mag_init(radius);mag_init(tail);arb_get_mag(radius,z.v);
        double upper=mag_get_d(radius);unsigned n=upper>0?unsigned(std::ceil((precision+32)/std::max(1.0,-std::log2(upper)))):2;
        n=std::max(2u,std::min(n,unsigned(precision+32)));Real term=z,sum;
        for(unsigned j=1;j<=n;++j){sum=sum+term/Real(long(j*(j+1)));term=-term*z;}
        mag_pow_ui(tail,radius,n+1);mag_mul_2exp_si(tail,tail,1);arb_add_error_mag(sum.v,tail);
        mag_clear(radius);mag_clear(tail);return sum;
    }
    return ((1+z)*lp(z)-z)/z;
}
Real triangle(const Real&u,const std::array<Real,3>&s){
    Real m=1+u*(1-u)*s[1],end=1+u*(1-u)*s[0],C=sq(1-u)*s[2],B=end-m+C;
    pos(m,"positive triangle start");pos(end,"positive triangle end");
    if(arb_is_zero(C.v))return real_log(m)+linear_mean(B/m);
    pos(C,"positive triangle quadratic coefficient");Real root=rt(sq(B)+4*C*m),alpha;
    if(arb_is_nonnegative(B.v))alpha=2*m/(root+B);
    else if(arb_is_negative(B.v))alpha=(root-B)/(2*C);
    else alpha=(root-B)/(2*C); // Both forms are exact; no sign-dependent mathematics.
    Real tau=2*end/(root-B+2*C);pos(alpha,"negative root outside simplex");pos(tau,"positive root outside simplex");
    return real_log(m)+linear_mean(1/alpha)+tau*lp(1/tau)-1;
}
Real value(const Graph&g,bool positive,const std::array<Real,6>&p){
    Real s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],u=p[5];
    Real t=sq(s/(1-s)),A=1+t,den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den);
    Real co,si;arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real radial=rt(t*z),sqrt3=rt(Real(3));
    std::array<std::array<int,3>,3>parent={g.parent[0],g.parent[1],g.parent[0]};
    for(int j=0;j<3;++j)parent[2][j]-=g.parent[1][j];
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)},ts0;
    Real V0;for(int j=0;j<3;++j){V0=V0+w[j]*Real(parent[j][0]*parent[j][0])*t;ts0[j]=Real(g.triangle[j][0]*g.triangle[j][0])*t;}
    Real base=(1-x)/x,at0=(1-x+V0)/x,a0=bubble_A(at0),K0=-(a0-bubble_A(base)),T0=-triangle(u,ts0),total;
    for(long sign:{-1L,1L}){
        Real qx=Real(sign)*radial*co,qy=Real(sign)*radial*si,q1=qx,q2=(-qx+sqrt3*qy)/2;
        auto norm=[&](const std::array<int,3>&r){
            if(positive)return sq(Real(r[0])*qx+Real(r[1])-Real(r[2])/2)+sq(Real(r[0])*qy+sqrt3*Real(r[2])/2)+Real(r[0]*r[0])*t*(1-z);
            return Real(r[0]*r[0])*t+2*Real(r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real V,lr;std::array<Real,3>ts;
        for(int j=0;j<3;++j){V=V+w[j]*norm(parent[j]);ts[j]=norm(g.triangle[j]);Real host=1+norm(g.host[j]);pos(host,"positive host propagator");lr=lr-lp((host-A)/A);}
        Real DK=-(bubble_A((1-x+V)/x)-a0),DT=-triangle(u,ts)-T0,ratio=real_exp(lr);
        total=total+(em1(lr)*K0*T0+ratio*(DK*T0+K0*DT+DK*DT))/2;
    }
    Real answer=2/(s*(1-s))*power3(t/A)*(1-z)*(1-x)*jac*(1-u)*total;
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite complete six-dimensional forest required");return answer;
}
#include "nested_closed_control_v1.hpp"
Real mapped(const Graph&g,bool positive,bool joint,bool control,const std::array<Real,6>&p){
    if(!joint)return controlled_value(g,positive,control,p);
    int index=g.parent[0][0]==g.parent[1][0]?0:1;
    if(index==1&&(g.parent[0][0]==0||g.parent[1][0]!=0))throw std::invalid_argument("unproved parent partition");
    for(int j=0;j<3;++j)if(g.triangle[j][0]*g.triangle[j][0]!=(j<2?1:0))throw std::invalid_argument("unproved triangle partition");
    Real t=sq(p[0]/(1-p[0])),ell=lp(t/2),xi=em1(p[3]*ell)/t,dx=ell*real_exp(p[3]*ell)/t;
    Real eta=em1(p[5]*ell)/t,du=ell*real_exp(p[5]*ell)/t,total;
    for(int side=0;side<2;++side){Real first=side?1-xi:xi;std::array<Real,3>bary;int at=0;bary[index]=first;
        for(int j=0;j<3;++j)if(j!=index){bary[j]=(1-first)*(at==0?p[4]:1-p[4]);++at;}
        Real x=bary[0],old=rt(x)/(rt(x)+rt(1-x)),den=sq(old)+sq(1-old),oldjac=2*old*(1-old)/sq(den);
        auto pp=p;pp[3]=old;pp[4]=bary[1]/(1-x);
        for(int end=0;end<2;++end){pp[5]=end?1-eta:eta;total=total+dx*(1-first)/((1-x)*oldjac)*du*controlled_value(g,positive,control,pp);}}
    return total;
}
Real radial_value(const Graph&g,bool positive,bool mixture,const std::array<Real,6>&p){
    Real s=p[0],x0=p[3],y=p[4],u=p[5],x=sq(x0)/(sq(x0)+sq(1-x0));
    std::array<Real,3>weights={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    std::array<int,3>loops={g.parent[0][0],g.parent[1][0],g.parent[0][0]-g.parent[1][0]};
    Real ell;for(int j=0;j<3;++j)ell=ell+weights[j]*Real(loops[j]*loops[j]);
    pos(ell,"positive parent radial scale");Real h=1/rt(rt(ell*u*(1-u))),den=1-s+h*s,w=h*s/den,jac=h/sq(den);
    auto pp=p;pp[0]=w;Real scaled=jac*controlled_value(g,positive,true,pp);
    if(!mixture)return scaled;
    Real ratio0=h/sq(h*(1-s)+s),ratio1=sq(den)/h;
    return controlled_value(g,positive,true,p)/(1+ratio0)+scaled*ratio1/(1+ratio1);
}
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original graph recipe");auto g=load_graph(id);if(g.sha!=sha)throw std::invalid_argument("original graph hash required");return g;}
}
int main(int argc,char**argv){try{using namespace nested_closed;if(argc!=2)throw std::invalid_argument("one recipe required");auto g=recipe(argv[1]);
    std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||(mode!="direct-original"&&mode!="direct-positive"&&mode!="joint-original"&&mode!="joint-positive"&&mode!="direct-control-original"&&mode!="direct-control-positive"&&mode!="joint-control-original"&&mode!="joint-control-positive"&&mode!="radial-control-original"&&mode!="radial-control-positive"&&mode!="radial-mixture-control-original"&&mode!="radial-mixture-control-positive")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded six-cube request required");
    precision=bits;bool positive=mode.find("positive")!=std::string::npos,joint=mode.rfind("joint-",0)==0,control=mode.find("control")!=std::string::npos;
    std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){std::array<Real,6>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open six-cube required");}if(i)std::cout<<',';interval(mode.rfind("radial-",0)==0?radial_value(g,positive,mode.find("mixture")!=std::string::npos,p):mapped(g,positive,joint,control,p));}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
