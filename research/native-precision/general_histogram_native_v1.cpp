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
#include "nested_five_functions_v1.hpp"
#include "g371_complete_functions_v1.hpp"
#include "general_two_host_functions_v1.hpp"
#include "general_two_host_angular_functions_v1.hpp"
Real general_angular_value(const DoubleHostGraph&g,bool positive,bool common_scale,const std::array<Real,7>&p){
    Real t,v,radial;
    if(common_scale){Real total=sq(p[0]/(1-p[0]));t=total*p[1];v=total*(1-p[1]);radial=sq(power3(total))*sq(p[1])*sq(1-p[1])/(2*p[0]*(1-p[0]));}
    else{t=sq(p[0]/(1-p[0]));v=sq(p[1]/(1-p[1]));radial=power3(t)*power3(v)/(p[0]*(1-p[0])*p[1]*(1-p[1]));}
    Real ct,st,cp,sp,cpsi,spsi;arb_sin_cos_pi(st.v,ct.v,(2*p[4]).v,precision);arb_sin_cos_pi(sp.v,cp.v,(2*p[5]).v,precision);arb_sin_cos_pi(spsi.v,cpsi.v,p[6].v,precision);
    Real z=p[2],Z=p[3];V4 e={rt(z)*ct,rt(z)*st,rt(1-z),Real(0)},n1={-st,ct,Real(0),Real(0)},n2={rt(1-z)*ct,rt(1-z)*st,-rt(z),Real(0)},n3={Real(0),Real(0),Real(0),Real(1)};
    V4 n=sum4(scale4(sum4(scale4(n1,cp),scale4(n2,sp)),rt(Z)),scale4(n3,rt(1-Z)));
    V4 q=scale4(e,rt(t)),r=scale4(sum4(scale4(e,cpsi),scale4(n,spsi)),rt(v));
    auto a=general_angular_chart(g,positive,q,r,2),b=general_angular_chart(g,positive,q,sum4(r,q,-1),1),c=general_angular_chart(g,positive,sum4(r,q,-1),q,0);
    Real residual=a[0]+b[0]+c[0],mean=a[1]+b[1]+c[1];Real orientation=3*(1-z)*rt(1-Z),intrinsic=Real(8)/3*sq(sq(spsi));
    Real answer=radial*intrinsic*(orientation*residual+mean);if(!arb_is_finite(answer.v))throw std::runtime_error("finite general angular control required");return answer;
}
#include "general_histogram_functions_v1.hpp"
}
int main(int argc,char**argv){try{using namespace nested_closed;if(argc!=2)throw std::invalid_argument("one recipe required");auto g=general_host_recipe(argv[1]);
    std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||(mode!="hist-pair-original"&&mode!="hist-pair-positive"&&mode!="hist-scale-original"&&mode!="hist-scale-positive"&&mode!="map-pair"&&mode!="map-scale")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded general seven-cube request required");
    precision=bits;bool positive=mode.find("positive")!=std::string::npos,scale=mode.find("scale")!=std::string::npos;
    if(mode=="map-pair"||mode=="map-scale"){
        std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"maps\":[";
        for(unsigned i=0;i<count;++i){std::array<Real,7>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict map test cube");}auto m=general_histogram_map(g.id,mode=="map-scale",p);if(i)std::cout<<',';std::cout<<"{\"point\":[";for(int j=0;j<7;++j){if(j)std::cout<<',';interval(m.point[j]);}std::cout<<"],\"jacobian\":";interval(m.jacobian);std::cout<<'}';}
        std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing map input");std::cout<<"]}\n";return 0;}
    std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){std::array<Real,7>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open seven-cube required");}if(i)std::cout<<',';interval(general_histogram_value(g,positive,scale,p));}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
