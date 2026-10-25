// Four unresolved E3 graphs: full paired forest, analytic common scale.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <algorithm>
#include "fmpq.h"
namespace weighted_e3 {
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
struct PrimitiveBase{Real d,t0,a0,b0;};
PrimitiveBase prepare(const Real&mass,const Real&d){pos(d,"positive weighted edge");pos(mass-d,"other mass positive");
 PrimitiveBase b;b.d=d;b.t0=(mass-d)/d;b.a0=bubble_A(b.t0);b.b0=bubble_B(b.t0);return b;}
std::array<Real,2> J_mass(const PrimitiveBase&b,const Real&z){
 if(arb_is_zero(z.v))return {Real(0),Real(0)};pos(z,"positive Schur momentum");Real t1=b.t0+z/b.d;
 return {b.a0-bubble_A(t1),(b.b0/sq(1+b.t0)-bubble_B(t1)/sq(1+t1))/b.d};
}
struct State{int sign;std::vector<unsigned> weighted;std::array<std::array<int,4>,6> routes;std::vector<std::array<int,2>> U,F;};
struct Graph{bool shared=false;std::string id,sha;std::vector<State> states;};
#include "weighted_e3_data.hpp"
Real det(const std::array<std::array<Real,3>,3>&a,unsigned n){
 std::array<unsigned,3> permutation={0,1,2};Real result;
 do{unsigned inversions=0;Real term(1);for(unsigned i=0;i<n;++i){term=term*a[i][permutation[i]];for(unsigned j=i+1;j<n;++j)inversions+=permutation[i]>permutation[j];}result=result+(inversions%2?-term:term);}
 while(std::next_permutation(permutation.begin(),permutation.begin()+n));return result;}
Real poly(const std::vector<std::array<int,2>>&terms,const std::array<Real,6>&alpha){
 Real total;for(const auto&term:terms){Real v(term[1]);for(unsigned j=0;j<6;++j)if(term[0]&(1<<j))v=v*alpha[j];total=total+v;}return total;}
std::array<Real,2> uf(const State&s,const std::array<Real,6>&alpha,bool positive){
 if(positive)return {poly(s.U,alpha),poly(s.F,alpha)};std::array<std::array<Real,4>,4>a;
 for(unsigned e=0;e<6;++e)for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)a[i][j]=a[i][j]+alpha[e]*Real(s.routes[e][i]*s.routes[e][j]);
 auto minor=[&](unsigned row,unsigned col){std::array<std::array<Real,3>,3>b;std::array<unsigned,3>ri={0,1,2+row},ci={0,1,2+col};
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)b[i][j]=a[ri[i]][ci[j]];return det(b,3);};
 return {a[0][0]*a[1][1]-a[0][1]*a[1][0],minor(0,0)+minor(1,1)-(minor(0,1)+minor(1,0))/2};}
Real distinct_density(const Graph&g,const std::array<Real,6>&alpha,const Real&tau,bool positive){
 Real mass,total;for(const auto&a:alpha)mass=mass+a;auto weights=g.states.at(0).weighted;if(weights.size()!=2)throw std::runtime_error("two distinct weights required");Real d=alpha[weights[0]],e=alpha[weights[1]];
 auto base=prepare(mass,d),shifted=prepare(mass+e*(1/tau-1),d);Real c=(tau-2)/sq(1-tau),h=1/(1-tau),k=1/(tau*sq(1-tau));
 for(const auto&s:g.states){if(s.weighted!=weights)throw std::runtime_error("fixed distinct weighted edges required");auto v=uf(s,alpha,positive);pos(v[0],"positive U");
   auto j=J_mass(base,v[1]/v[0]),other=J_mass(shifted,v[1]/v[0]);total=total+Real(s.sign)*(c*j[0]-h*e*j[1]+k*other[0])/power3(v[0]);}
 return total;}
Real shared_density(const Graph&g,const std::array<Real,6>&alpha,const Real&v,bool positive){
 Real mass,total;for(const auto&a:alpha)mass=mass+a;auto weights=g.states.at(0).weighted;
 if(weights.size()!=1)throw std::runtime_error("one shared weight required");Real d=alpha[weights[0]];
 for(const auto&s:g.states){if(s.weighted!=weights)throw std::runtime_error("fixed shared edge required");auto pol=uf(s,alpha,positive);pos(pol[0],"positive U");
   Real z=pol[1]/pol[0],a=mass+v*z,b=bubble_B((a-d)/d);total=total+-Real(s.sign)*sq(d)*z*sq(b)/(power3(a)*power3(pol[0]));}
 return total;}
Real value(const Graph&g,const std::array<Real,6>&p,bool positive,bool square){
 std::array<Real,5>x;Real jac(1),total;for(unsigned i=0;i<5;++i){x[i]=square?sq(p[i]):p[i];if(square)jac=jac*2*p[i];}
 Real r=1-sq(p[5]),tau=(1-sq(r))/4,spectral=g.shared?Real(1):p[5]*sq(sq(r))/6;
 for(unsigned pivot=0;pivot<6;++pivot){std::array<Real,6>alpha;unsigned next=0;for(unsigned e=0;e<6;++e)alpha[e]=e==pivot?Real(1):x[next++];
  total=total+(g.shared?shared_density(g,alpha,p[5],positive):distinct_density(g,alpha,tau,positive));}
 total=total*jac*spectral;if(!arb_is_finite(total.v))throw std::runtime_error("finite full E3 forest required");return total;}
#include "weighted_e3_mobius_functions_v1.hpp"
#include "shared_bubble_square_native_v1.hpp"
#include "shared_e3_closed_functions_v1.hpp"
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one double weighted E2 recipe required");
 auto g=load_weighted_e2(id);if(g.sha!=sha)throw std::invalid_argument("original double weighted identity required");return g;}
}
int main(int argc,char**argv){try{using namespace weighted_e3;
 if(argc!=2)throw std::invalid_argument("one E3 recipe required");auto graph=recipe(argv[1]);std::string mode;long bits,count;
 if(!(std::cin>>mode>>bits>>count)||(mode!="sector2-original"&&mode!="sector2-positive"&&mode!="sector4-original"&&mode!="sector4-positive"&&mode!="map-sector2"&&mode!="map-sector4")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded exact six-dimensional request required");
 bool fourth=mode.find("sector4")!=std::string::npos,positive=mode.find("positive")!=std::string::npos,maps=mode.find("map-")==0;
 std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\""<<(maps?"maps":"values")<<"\":[";
 for(long i=0;i<count;++i){std::array<Rat,6>point;for(auto&v:point){v=read_rational();if(!(Rat(0)<v)||!(v<Rat(1)))throw std::invalid_argument("strict open exact six-cube required");}
  if(i)std::cout<<',';precision=bits;
  if(maps){std::array<Rat,6>sector;for(long j=0;j<6;++j)sector[j]=point[j];rational_sector_map_output(e3_mobius_table(graph.id,fourth),sector);}
  else interval(shared_closed_value(graph,point,positive,fourth,bits));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing E3 sector input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
