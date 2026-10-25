// Frozen graph identity and two independent exact auxiliary representations.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace nested_bubble_vertex {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
namespace M=disjoint_v_integrated;
Real sq(const Real&x){Real r;arb_sqr(r.v,x.v,precision);return r;}
Real rt(const Real&x){Real r;arb_sqrt(r.v,x.v,precision);return r;}
Real lp(const Real&x){Real r;arb_log1p(r.v,x.v,precision);return r;}
Real ex(const Real&x){Real r;arb_exp(r.v,x.v,precision);return r;}
Real em1(const Real&x){Real r;arb_expm1(r.v,x.v,precision);return r;}
Real ah(const Real&x){Real r;arb_atanh(r.v,x.v,precision);return r;}
void pos(const Real&x,const char*why){if(!arb_is_positive(x.v))throw std::runtime_error(why);}
Real bubble(const Real&t,bool primitive){
    pos(t,"positive bubble argument");
    if(primitive){auto b=disjoint_interval::bubble(fourier_e2_root::Z(t));
        if(!arb_contains_zero(acb_imagref(b.v)))throw std::runtime_error("real bubble required");
        Real out;arb_set(out.v,acb_realref(b.v));return out;}
    Real z=rt(t/(t+4));pos(1-z,"independent bubble atanh domain");
    return (t+4)/3*(ah(z)/z-1)-t/9;
}
Real bubble_prime(const Real&t){
    pos(t,"positive bubble derivative argument");Real z=rt(t/(t+4));pos(1-z,"bubble derivative atanh domain");
    return ((t-2)*ah(z)/z+2)/(3*t)-Real(5)/18;
}
Real value(const std::string&mode,bool control,const std::array<Real,6>&p){
    const auto&s=p[0];const auto&z=p[1];const auto&theta=p[2];const auto&x0=p[3];const auto&y=p[4];const auto&u=p[5];
    bool primitive=mode=="primitive";
    Real t=sq(s/(1-s)),A=1+t,X=2+t,den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den);
    Real co,si;arb_sin_cos_pi(si.v,co.v,theta.v,precision);
    Real p2=(-co+rt(Real(3))*si)/2,q=rt(t*z),l1=-2*q*co,l2=2*q*p2;
    Real alpha=x*(1-x),V0=alpha*t,L=2*alpha*q*(y*p2-(1-y)*co),C=alpha+sq(1-x)*y*(1-y);
    Real K0,a,d,base,w;
    if(primitive){Real ell=lp(V0),arg=(1-x+em1(u*ell))/x;K0=-ell*bubble(arg,true)/(1+arg);}
    else{
        Real b=1-x+V0,ell=lp(b/(4*x));a=x/b*em1(u*ell);d=1-a;base=x+a*(1-x);
        pos(a,"positive spectral parameter");pos(1-4*a,"spectral interval interior");
        w=(1-4*a)*rt(1-4*a)/6*x/b*ell*ex(u*ell);
        K0=w*((1+d)/sq(d)*lp(V0)-x/d*V0/(1+V0)-lp(a*V0/base)/(a*sq(d)));
    }
    Real total;
    for(long sign:{-1L,1L}){
        Real DV=C+sign*L,V=V0+DV;pos(V,"positive child variance");
        Real denomin=(X+sign*l1)/A*(X+sign*l2)/A;pos(denomin,"positive host product");
        Real change=-(2/A+1/sq(A)+sign*(X/A)*(l1+l2)/A+(l1/A)*(l2/A))/denomin,delta;
        if(primitive){Real shift=lp(DV/(1+V0)),arg=((1+V0)*ex(u*shift)-x)/x;
            delta=-shift*bubble(arg,true)/(1+arg);
        }else delta=w*((1+d)/sq(d)*lp(DV/(1+V0))-x/d*DV/((1+V0)*(1+V))-lp(a*DV/(base+a*V0))/(a*sq(d)));
        total=total+(change*K0+delta/denomin)/2;
    }
    Real common=2/(s*(1-s))*M::power(t/A,3)*(1-x)*jac;
    Real result=common*(1-z)*total;
    if(control){
        Real q1=q*co,q2=q*p2,arg=(1-x+V0)/x,B=bubble(arg,primitive);
        Real first=-x*B/sq(1+V0),second=2*x*B/M::power(1+V0,3)-bubble_prime(arg)/sq(1+V0);
        Real point=(4*(sq(q1)+sq(q2)-q1*q2)/sq(A)-2/A)*K0+C*first+sq(L)/2*second-2*(q2-q1)*L/A*first;
        Real average=(5*t/(3*sq(A))-2/A)*K0+(alpha+sq(1-x)/6-t*alpha/A)*first+5*t*sq(alpha)/18*second;
        result=result+common*(average/2-(1-z)*point);
    }
    if(!arb_is_finite(result.v))throw std::runtime_error("finite nested full kernel required");return result;
}
void recipe(const char*path){
    std::ifstream in(path);std::string id,sha,extra;
    if(!(in>>id>>sha)||id!="L3E3G8"||sha!="6032d8167eb696df1d3559d8d2de0920cec0012d19b90a6d2fab8cdaccbc5ee5"||(in>>extra))
        throw std::invalid_argument("frozen nested graph recipe required");
}
}
int main(int argc,char**argv){try{using namespace nested_bubble_vertex;
    if(argc!=2)throw std::invalid_argument("one frozen recipe required");recipe(argv[1]);
    std::string mode;unsigned bits,orbit,control,count;
    if(!(std::cin>>mode>>bits>>orbit>>control>>count)||(mode!="primitive"&&mode!="spectral")||bits<512||bits>2048||
       (orbit!=1&&orbit!=4&&orbit!=8)||control>1||count<1||count>512)throw std::invalid_argument("bounded nested batch required");
    precision=bits;
    std::cout<<"{\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"orbit\":"<<orbit<<",\"control\":"<<control<<",\"auxiliary_antithetic\":2,\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){
        std::array<Real,6>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open six-cube required");}
        Real sum;
        for(unsigned j=0;j<orbit;++j){auto pp=p;pp[2]=p[2]+Real(long(j))/Real(long(orbit));
            if(arb_ge(pp[2].v,Real(1).v))pp[2]=pp[2]-1;
            else if(!arb_lt(pp[2].v,Real(1).v))throw std::runtime_error("uncertain angular orbit boundary");
            Real one=value(mode,control,pp);pp[5]=1-pp[5];
            sum=sum+(one+value(mode,control,pp))/(2*Real(long(orbit)));}
        if(i)std::cout<<',';interval(sum);
    }
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing nested request");
    std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
