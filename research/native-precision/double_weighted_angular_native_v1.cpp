// Complete G9 angular integrals: closed logarithm/dilogarithm and bounded positive series.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace double_angular {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
#include "shared_overlap_identity.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real cube(const Real&a){return a*sq(a);}
struct Shape{Real k,r;};
Shape shape(const Real&s,const Real&t){Real a=1+s+t,disc=rt(sq(a)-4*s*t),k=2/(a+disc);return {k,k*rt(s*t)};}
Real f0(const Shape&a){return a.k*(1-sq(a.r)/3);}
Real f1(const Shape&a){return a.k*a.r*(Real(1)/2-sq(a.r)/4);}
Real angular(const Real&s,const Real&t,bool series){
 Shape a=shape(s,Real(1)),b=shape(t,Real(1)),c=shape(s,t);Real z=a.r*b.r*c.r;
 if(!arb_is_positive(z.v)||!arb_lt(z.v,(Real(7)/40).v))throw std::runtime_error("positive uniformly bounded angular series required");
 Real aa=sq(a.r),bb=sq(b.r),cc=sq(c.r),prefactor=a.k*b.k*c.k;
 if(series){
  constexpr long N=256;Real sum(0),zn(1);
  for(long n=0;n<N;++n){
   Real factor=Real(2*(n+1)*(n+3))/Real(3*(n+2));
   sum=sum+zn*factor*(Real(1)/Real(n+1)-aa/Real(n+3))*(Real(1)/Real(n+1)-bb/Real(n+3))*(Real(1)/Real(n+1)-cc/Real(n+3));zn=zn*z;
  }
  Real bound=zn/(Real((N+1)*(N+1))*(1-z)),tail,zero(0);arb_union(tail.v,zero.v,bound.v,precision);
  return prefactor*(sum+tail);
 }
 Real e1=aa+bb+cc,e2=aa*bb+aa*cc+bb*cc,e3=aa*bb*cc,log,dilog;
 arb_log1p(log.v,(-z).v,precision);log=-log;arb_polylog_si(dilog.v,2,z.v,precision);
 Real one1=log/z,two1=dilog/z,one2=(log-z)/sq(z),one3=(log-z-sq(z)/2)/cube(z),two3=(dilog-z-sq(z)/4)/cube(z);
 return prefactor*Real(2)/3*(2*two1-(1+e1)*one1+(1+e1+e2+e3)*one2-(e2+e3)*one3-2*e3*two3);
}
Real bubble(const Real&s){Real h;arb_asinh(h.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*h/3-4*s/9-Real(4)/3;}
Real remainder(const Real&s,const Real&t,bool series){
 Real ds=1+s,dt=1+t;Shape h=shape(s,t);
 Real js=-Real(1)/sq(ds)+2*s/(3*cube(ds)),jt=-Real(1)/sq(dt)+2*t/(3*cube(dt));
 Real mass=f0(h)/(ds*dt),kinetic=f0(h)*(js/dt+jt/ds)+4*rt(s*t)*f1(h)/(9*sq(ds)*sq(dt));
 Real ps=(f0(shape(s,Real(1)))-1/ds-js)/sq(dt),pt=(f0(shape(t,Real(1)))-1/dt-jt)/sq(ds);
 return angular(s,t,series)-mass-kinetic-ps-pt;
}
Real density(const std::array<Real,2>&p,bool series,bool cubic){
 Real dx=1-p[0],dy=1-p[1],s=p[0]/(cubic?cube(dx):dx),t=p[1]/(cubic?cube(dy):dy);
 Real ds=(cubic?1+2*p[0]:Real(1))/(cubic?sq(sq(dx)):sq(dx));
 Real dt=(cubic?1+2*p[1]:Real(1))/(cubic?sq(sq(dy)):sq(dy));
 Real value=ds*dt*sq(s)*sq(t)*bubble(s)*bubble(t)/(4*sq(1+s)*sq(1+t))*remainder(s,t,series);
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete angular G9 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G9 recipe required");}
}
int main(int argc,char**argv){try{using namespace double_angular;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;
 if(!(std::cin>>mode>>bits>>count)||
 (mode!="angular-two-linear-closed"&&mode!="angular-two-linear-series"&&mode!="angular-two-cubic-closed"&&mode!="angular-two-cubic-series")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded two-dimensional angular request required");
 precision=bits;bool series=mode.find("series")!=std::string::npos,cubic=mode.find("cubic")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,2>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open two-cube required");}
 if(i)std::cout<<',';interval(density(p,series,cubic));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing angular request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
