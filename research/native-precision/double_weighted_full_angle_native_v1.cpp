// G20 full external angular average; bounded Gegenbauer and Chebyshev recurrences.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace full_angle {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
#include "shared_overlap_identity.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real absolute(const Real&a){Real r;arb_abs(r.v,a.v);return r;}
Real cube(const Real&a){return a*sq(a);}
Real symmetric(const Real&a){Real r;arb_union(r.v,(-a).v,a.v,precision);return r;}
struct Shape{Real k,r,g;};
Shape shape(const Real&s,const Real&e){
 Real a=1+s+e,delta=rt(sq(a)-4*s*e),k=2/(a+delta);
 return {k,k*rt(s*e),-(1+(a-2*s)/delta)/(a+delta)};
}
Real f0(const Shape&a){return a.k*(1-sq(a.r)/3);}
Real pair(const Real&s,const Real&t,const Real&dot,const Real&e,bool derivative,bool chebyshev){
 Shape a=shape(s,e),b=shape(t,e);Real pref=a.k*b.k,z=a.r*b.r;
 if(!arb_is_positive(z.v)||!arb_lt(z.v,(Real(1)/4).v))throw std::runtime_error("uniform pair tail ratio required");
 Real q=pref*e*dot,z2=sq(pref*e)*s*t,aa=sq(a.r),bb=sq(b.r),g=a.g+b.g;
 Real old2(1),old1=4*q,uold2(1),uold1=2*q,total(0),zn(1);
 constexpr long N=256;
 for(long n=0;n<N;++n){
  Real poly;
  if(n==0)poly=old2;
  else if(n==1)poly=old1;
  else if(chebyshev){Real unew=2*q*uold1-z2*uold2;poly=Real(n+1)*unew+z2*old2;uold2=uold1;uold1=unew;old2=old1;old1=poly;}
  else {poly=(Real(2*(n+1))*q*old1-Real(n+2)*z2*old2)/Real(n);old2=old1;old1=poly;}
  Real da=Real(1)/Real(n+1)-aa/Real(n+3),db=Real(1)/Real(n+1)-bb/Real(n+3),factor=da*db;
  if(derivative){
   Real dap=-2*aa*(a.g+1/(2*e))/Real(n+3),dbp=-2*bb*(b.g+1/(2*e))/Real(n+3);
   factor=(Real(n+1)*g+Real(n)/e)*da*db+dap*db+da*dbp;
  }
  total=total+2*poly*factor/Real(n+2);zn=zn*z;
 }
 Real bound=derivative?zn*(Real(4*N+10)/(1-z)+4*z/sq(1-z)):zn/(1-z);
 return pref*(total+symmetric(bound));
}
Real pair_second_bound(const Real&s,const Real&t){
 Real a=rt(2*s),b=rt(2*t);return 2*sq(1+a)+2*sq(1+b)+2*(1+a)*(1+b)+a+b;
}
Real angular(const Real&s,const Real&t,const Real&h,bool chebyshev){
 Real third=s+t-2*h,c=2+2*h,e(1);bool derivative=arb_contains_zero(c.v);
 if(derivative){
  e=1-c;
  if(!arb_gt(e.v,(Real(1)/2).v)||!arb_lt(e.v,(Real(3)/2).v))throw std::runtime_error("certified removable-limit range required");
 }
 Real value=pair(s,third,s-h,e,derivative,chebyshev)+pair(t,third,t-h,e,derivative,chebyshev)-pair(s,t,-h,e,derivative,chebyshev);
 if(!derivative)return value/c;
 Real bound=absolute(c)*(pair_second_bound(s,third)+pair_second_bound(t,third)+pair_second_bound(s,t))/2;
 std::cerr<<"full-angle-removable-limit\n";
 return value+symmetric(bound);
}
Real bubble(const Real&s){Real a;arb_asinh(a.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*a/3-4*s/9-Real(4)/3;}
Real remainder(const Real&s,const Real&t,const Real&h,bool chebyshev){
 Real third=s+t-2*h,ds=1+s,dt=1+t,dl=1+third;
 Real mass=1/(ds*dt*dl);
 Real kinetic=mass*(-1/ds-1/dt-1/dl+Real(2)/3*(s/sq(ds)+t/sq(dt)+third/sq(dl)-h/(ds*dt)+(s-h)/(ds*dl)+(t-h)/(dt*dl)));
 Real js=-Real(1)/sq(ds)+2*s/(3*cube(ds)),jt=-Real(1)/sq(dt)+2*t/(3*cube(dt));
 Real ps=(f0(shape(s,Real(1)))-1/ds-js)/sq(dt),pt=(f0(shape(t,Real(1)))-1/dt-jt)/sq(ds);
 return angular(s,t,h,chebyshev)-mass-kinetic-ps-pt;
}
Real density(const std::array<Real,3>&p,bool chebyshev,bool cubic){
 Real dx=1-p[0],dy=1-p[1],s=p[0]/(cubic?cube(dx):dx),t=p[1]/(cubic?cube(dy):dy);
 Real ds=(cubic?1+2*p[0]:Real(1))/(cubic?sq(sq(dx)):sq(dx));
 Real dt=(cubic?1+2*p[1]:Real(1))/(cubic?sq(sq(dy)):sq(dy));
 Real pi,co,si;arb_const_pi(pi.v,precision);arb_cos_pi(co.v,p[2].v,precision);arb_sin_pi(si.v,p[2].v,precision);
 Real h=rt(s*t)*co,weight=Real(8)/3*sq(sq(si));
 Real value=ds*dt*sq(s)*sq(t)*bubble(s)*bubble(t)/(4*sq(1+s)*sq(1+t))*weight*remainder(s,t,h,chebyshev);
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete full-angle G20 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G20 recipe required");}
}
int main(int argc,char**argv){try{using namespace full_angle;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;
 if(!(std::cin>>mode>>bits>>count)||
 (mode!="full-angle-linear-gegenbauer"&&mode!="full-angle-linear-chebyshev"&&mode!="full-angle-cubic-gegenbauer"&&mode!="full-angle-cubic-chebyshev")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded three-dimensional full-angle request required");
 precision=bits;bool chebyshev=mode.find("chebyshev")!=std::string::npos,cubic=mode.find("cubic")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,3>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open three-cube required");}
 if(i)std::cout<<',';interval(density(p,chebyshev,cubic));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing full-angle request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
