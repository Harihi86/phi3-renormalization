// Twenty-two original one-bubble E3 graphs, complete vertex forests and common scale.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <algorithm>
#include "fmpq.h"
namespace single_weighted_e3 {
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::Z;
using fourier_e2_root::readq;using fourier_e2_root::interval;
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

Real weighted_C(const Real&t){
    if(arb_is_zero(t.v))return Real(0);pos(t,"nonnegative weighted primitive argument");
    Real h;arb_asinh(h.v,(rt(t)/2).v,precision);
    return h*rt(t*(t+4))*(t-2)/(3*(t+1))+4*sq(h)/3-Real(11)*t/18+lp(t)/18+8*t/(9*(t+1));
}

Real bubble_B(const Real&t){
 if(arb_is_zero(t.v))return Real(0);pos(t,"positive bubble argument");Real h;arb_asinh(h.v,(rt(t)/2).v,precision);
 return (t+4)*rt(t+4)/rt(t)*h/3-4*t/9-Real(4)/3;
}

struct State{int sign;unsigned weighted;std::array<std::array<int,5>,9>routes;std::vector<std::array<int,2>>U,F;};
struct Graph{std::string id,sha;std::vector<State>states;};
#include "single_weighted_e3_data.hpp"
#include "single_weighted_e3_functions_v1.hpp"
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one single-bubble E3 recipe required");
 auto g=load_single_weighted_e3(id);if(g.sha!=sha)throw std::invalid_argument("original graph identity required");return g;}
}
int main(int argc,char**argv){try{using namespace single_weighted_e3;
 if(argc!=2)throw std::invalid_argument("one original E3 recipe required");auto graph=recipe(argv[1]);std::string mode;long bits,count;
 if(!(std::cin>>mode>>bits>>count)||(mode!="sector2-original"&&mode!="sector2-positive"&&mode!="sector4-original"&&mode!="sector4-positive"&&mode!="map-sector2"&&mode!="map-sector4")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded exact nine-dimensional request required");
 bool fourth=mode.find("sector4")!=std::string::npos,positive=mode.find("positive")!=std::string::npos,maps=mode.find("map-")==0;
 std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\""<<(maps?"maps":"values")<<"\":[";
 for(long i=0;i<count;++i){std::array<Rat,9>point;for(auto&v:point){v=read_rational();if(!(Rat(0)<v)||!(v<Rat(1)))throw std::invalid_argument("strict open exact nine-cube required");}
  if(i)std::cout<<',';precision=bits;
  if(maps)rational_sector_map_output(single_sector_table(graph.id,fourth),point);
  else interval(single_value(graph,point,positive,fourth,bits));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing nine-cube input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
