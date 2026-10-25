// Original G8 forests with full spherical pair and shifted weighted-edge Taylor jets.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace weighted_pair {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
#include "shared_overlap_identity.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real ab(const Real&a){Real r;arb_abs(r.v,a.v);return r;}
Real cube(const Real&a){return a*sq(a);}
Real sym(const Real&a){Real r;arb_union(r.v,(-a).v,a.v,precision);return r;}
struct Shape{Real k,r;};
Shape shape(const Real&s,const Real&e){Real a=1+s+e,k=2/(a+rt(sq(a)-4*s*e));return {k,k*rt(s*e)};}
Real hfunction(const Real&z,bool asin_kernel){
 if(arb_lt(ab(z).v,(Real(1)/4).v)){
  Real total(0),power(1),tolerance(1);arb_mul_2exp_si(tolerance.v,tolerance.v,-long(precision));
  for(long n=0;n<4096;++n){
   total=total+power/Real((2*n+1)*(2*n+3));power=-power*z;
   Real tail=ab(power)/(Real((2*n+3)*(2*n+5))*(1-ab(z)));
   if(arb_lt(tail.v,tolerance.v))return total+sym(tail);
  }
  throw std::runtime_error("bounded H-series termination required");
 }
 if(!arb_is_positive(z.v))throw std::runtime_error("positive H argument required");
 Real angle,root=rt(z);
 if(asin_kernel)arb_asin(angle.v,rt(z/(1+z)).v,precision);
 else arb_atan(angle.v,root.v,precision);
 return ((1+z)*angle/root-1)/(2*z);
}
Real pair(const Real&t,const Real&u,const Real&dot,const Real&gram,const Real&s,bool asin_kernel){
 Shape a=shape(t,s),b=shape(u,s);Real pref=a.k*b.k,q=pref*s*dot,eta2=sq(pref*s)*gram;
 Real denominator=1-q,v=sq(a.r)+sq(b.r)-2*q;
 return pref*(1-v/denominator*hfunction(eta2/sq(denominator),asin_kernel));
}
Real f0(const Shape&a){return a.k*(1-sq(a.r)/3);}
Real bubble(const Real&s){Real a;arb_asinh(a.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*a/3-4*s/9-Real(4)/3;}
Real subtractions(const Real&s,const Real&t){
 Real ds=1+s,dt=1+t;Shape h=shape(s,t);Real k=h.k,r2=sq(h.r);
 Real j1=k*(1-r2/3),j2=sq(k),j3=cube(k)/(1-r2);
 Real a=Real(2)/3*(-1/ds-1/dt+s/sq(ds)+t/sq(dt))-(1+s+t)/(3*ds*dt);
 Real b=-Real(1)/3+(s-t-1)/(3*ds)+(t-s-1)/(3*dt);
 Real mass=j1/(ds*dt),kinetic=(a*j1+b*j2-2*j3/3+1/(3*ds*dt))/(ds*dt);
 Real js=-Real(1)/sq(ds)+2*s/(3*cube(ds)),jt=-Real(1)/sq(dt)+2*t/(3*cube(dt));
 Real ps=(f0(shape(s,Real(1)))-1/ds-js)/sq(dt),pt=(f0(shape(t,Real(1)))-1/dt-jt)/sq(ds);
 return mass+kinetic+ps+pt;
}

Real bare(const Real&s,const Real&t,const Real&b,bool asin_kernel){
 Real root=rt(t),u=1+t+2*root*b,dot=1+root*b;
 return pair(Real(1),u,dot,t*(1-sq(b)),s,asin_kernel)/(1+u);
}
Real angular(const Real&s,const Real&t,bool asin_kernel){
 Real sub=subtractions(s,t),strip_s=2+s-Real(5)/2*rt(s),strip_t=2+t-Real(5)/2*rt(t);
 for(long n=32;n<=2048;n*=2){
  Real total(0);
  for(long j=0;j<n;++j){
   Real fraction=Real(2*j+1)/Real(2*n),b,si;arb_cos_pi(b.v,fraction.v,precision);arb_sin_pi(si.v,fraction.v,precision);
   total=total+Real(8)/3*sq(sq(si))*bare(s,t,b,asin_kernel);
  }
  Real value=total/Real(n)-sub,power4(1);arb_mul_2exp_si(power4.v,power4.v,2*n);
  Real error=Real(1250)/(42*strip_s*strip_t*(power4-1)),tolerance=ab(value);
  arb_mul_2exp_si(tolerance.v,tolerance.v,-80);
  if(arb_lt(error.v,tolerance.v)){
   std::cerr<<"cauchy-angle-count "<<n<<'\n';return value+sym(error);
  }
 }
 throw std::runtime_error("complete angular Cauchy enclosure did not meet fixed relative tolerance");
}
Real density(const std::array<Real,2>&p,bool asin_kernel,bool diagonal){
 Real dx=1-p[0],dy=1-p[1],a=p[0]/cube(dx),b=p[1]/cube(dy),da=(1+2*p[0])/sq(sq(dx)),db=(1+2*p[1])/sq(sq(dy)),value;
 if(!diagonal)value=da*db*sq(a)*sq(b)*bubble(a)*bubble(b)/(4*sq(1+a)*sq(1+b))*angular(a,b,asin_kernel);
 else{Real s=sq(a+b),t=sq(a),jac=4*a*(a+b)*da*db;value=jac*sq(s)*sq(t)*bubble(s)*bubble(t)/(4*sq(1+s)*sq(1+t))*(angular(s,t,asin_kernel)+angular(t,s,asin_kernel));}
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete original G20 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G20 recipe required");}
}
int main(int argc,char**argv){try{using namespace weighted_pair;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);std::string mode;unsigned bits,count;
 if(!(std::cin>>mode>>bits>>count)||(mode!="cauchy-cubic-atan"&&mode!="cauchy-cubic-asin"&&mode!="cauchy-diagonal-atan"&&mode!="cauchy-diagonal-asin")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded complete G20 request required");
 precision=bits;bool asin_kernel=mode.find("asin")!=std::string::npos,diagonal=mode.find("diagonal")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,2>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open two-cube required");}
 if(i)std::cout<<',';interval(density(p,asin_kernel,diagonal));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing G20 request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
