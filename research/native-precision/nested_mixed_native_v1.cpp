// Complete eight-forest mixed graphs after three analytic child-loop integrations.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace nested_mixed {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
namespace M=disjoint_v_integrated;
struct Graph{std::string id,sha;std::array<std::array<int,3>,2> parent;std::array<std::array<int,3>,3> triangle,host;};
#include "nested_mixed_data.hpp"
Real sq(const Real&x){Real r;arb_sqr(r.v,x.v,precision);return r;}
Real rt(const Real&x){Real r;arb_sqrt(r.v,x.v,precision);return r;}
Real lp(const Real&x){Real r;arb_log1p(r.v,x.v,precision);return r;}
Real ex(const Real&x){Real r;arb_exp(r.v,x.v,precision);return r;}
Real em1(const Real&x){Real r;arb_expm1(r.v,x.v,precision);return r;}
Real ah(const Real&x){Real r;arb_atanh(r.v,x.v,precision);return r;}
void pos(const Real&x,const char*why){if(!arb_is_positive(x.v))throw std::runtime_error(why);}
Real bubble(const Real&t,bool primitive){
    pos(t,"positive subtracted bubble argument");
    if(primitive){auto b=disjoint_interval::bubble(fourier_e2_root::Z(t));
        if(!arb_contains_zero(acb_imagref(b.v)))throw std::runtime_error("real subtracted bubble required");
        Real out;arb_set(out.v,acb_realref(b.v));return out;}
    Real z=rt(t/(t+4));pos(1-z,"independent bubble atanh domain");
    return (t+4)/3*(ah(z)/z-1)-t/9;
}
std::pair<Real,Real> local(bool primitive,const Real&x,const Real&u,const Real&V0,const Real&DV){
    if(primitive){Real ell=lp(V0),a=(1-x+em1(u*ell))/x;
        Real K0=-ell*bubble(a,true)/(1+a);
        Real shift=lp(DV/(1+V0)),b=((1+V0)*ex(u*shift)-x)/x;
        return {K0,-shift*bubble(b,true)/(1+b)};}
    Real b=1-x+V0,ell=lp(b/(4*x)),a=x/b*em1(u*ell),d=1-a,base=x+a*(1-x);
    pos(a,"positive spectral coordinate");pos(1-4*a,"interior bubble spectrum");
    Real weight=(1-4*a)*rt(1-4*a)/6*x/b*ell*ex(u*ell);
    Real K0=weight*((1+d)/sq(d)*lp(V0)-x/d*V0/(1+V0)-lp(a*V0/base)/(a*sq(d)));
    Real DK=weight*((1+d)/sq(d)*lp(DV/(1+V0))-x/d*DV/((1+V0)*(1+V0+DV))-lp(a*DV/(base+a*V0))/(a*sq(d)));
    return {K0,DK};
}
Real value(const Graph&g,const std::string&mode,const std::array<Real,8>&p){
    Real s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],u=p[5];bool primitive=mode=="primitive";
    Real t=sq(s/(1-s)),A=1+t,den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den);
    Real co,si;arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real radial=rt(t*z);
    Real rem=rt(1-p[6]),v0=p[6]/(1+rem),v1=rem*p[7],v2=rem*(1-p[7]);
    std::array<Real,3>pw={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)},tw={v0*v1,v0*v2,v1*v2};
    std::array<std::array<int,3>,3>parent={g.parent[0],g.parent[1],g.parent[0]};
    for(int j=0;j<3;++j)parent[2][j]-=g.parent[1][j];
    Real V0,VT0;for(int j=0;j<3;++j){V0=V0+pw[j]*Real(parent[j][0]*parent[j][0])*t;VT0=VT0+tw[j]*Real(g.triangle[j][0]*g.triangle[j][0])*t;}
    Real total;
    for(long sign:{-1L,1L}){
        Real q1=Real(sign)*radial*co,q2=Real(sign)*radial*(-co+rt(Real(3))*si)/2;
        auto change=[&](const std::array<int,3>&r){return Real(2*r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real DV,DVT,logratio;
        for(int j=0;j<3;++j){DV=DV+pw[j]*change(parent[j]);DVT=DVT+tw[j]*change(g.triangle[j]);
            Real shift=change(g.host[j])/A;pos(1+shift,"positive physical host denominator");logratio=logratio-lp(shift);}
        pos(1+V0+DV,"positive nested vertex mass");pos(1+VT0+DVT,"positive disjoint triangle mass");
        auto nested=local(primitive,x,u,V0,DV);Real K0=nested.first,DK=nested.second,T0=-lp(VT0),DT=-lp(DVT/(1+VT0));
        total=total+(em1(logratio)*K0*T0+ex(logratio)*(DK*T0+K0*DT+DK*DT))/2;
    }
    // Triangle equal-area map contributes 1/2; the parent uses its own Jacobian.
    Real result=2/(s*(1-s))*M::power(t/A,3)*(1-z)*(1-x)*jac/2*total;
    if(!arb_is_finite(result.v))throw std::runtime_error("finite full mixed forest required");return result;
}
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
    if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one frozen mixed graph recipe required");
    auto g=load_graph(id);if(g.sha!=sha)throw std::invalid_argument("original graph digest required");return g;}
}
int main(int argc,char**argv){try{using namespace nested_mixed;
    if(argc!=2)throw std::invalid_argument("one graph recipe required");auto g=recipe(argv[1]);
    std::string mode;unsigned bits,orbit,control,count;
    if(!(std::cin>>mode>>bits>>orbit>>control>>count)||(mode!="primitive"&&mode!="spectral")||bits<512||bits>2048||
       (orbit!=1&&orbit!=4&&orbit!=8)||control!=0||count<1||count>512)throw std::invalid_argument("bounded complete mixed batch required");
    precision=bits;
    std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"orbit\":"<<orbit<<",\"control\":"<<control<<",\"auxiliary_antithetic\":2,\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){std::array<Real,8>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open eight-cube required");}
        Real sum;for(unsigned j=0;j<orbit;++j){auto pp=p;pp[2]=p[2]+Real(long(j))/Real(long(orbit));
            if(arb_ge(pp[2].v,Real(1).v))pp[2]=pp[2]-1;
            else if(!arb_lt(pp[2].v,Real(1).v))throw std::runtime_error("uncertain exact angular orbit boundary");
            Real a=value(g,mode,pp);pp[5]=1-pp[5];sum=sum+(a+value(g,mode,pp))/(2*Real(long(orbit)));}
        if(i)std::cout<<',';interval(sum);}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing mixed request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
