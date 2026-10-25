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
Real bubble(const Real&s){Real a;arb_asinh(a.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*a/3-4*s/9-Real(4)/3;}
Real weighted(const Real&s){return bubble(s)/sq(1+s);}
struct Weight{Real w,first,second;};
Weight weight_jet(const Real&s){
 Real d=1+s,a;arb_asinh(a.v,(rt(s)/2).v,precision);Real b=bubble(s),j=4*a/rt(s*(s+4));
 Real bp=rt(s+4)*(s-2)*a/(3*s*rt(s))+(s+4)/(6*s)-Real(4)/9,bpp=1/(6*s)-(1-j)/sq(s);
 return {b/sq(d),bp/sq(d)-2*b/cube(d),bpp/sq(d)-4*bp/cube(d)+6*b/sq(sq(d))};
}
Real subtractions(const Real&s,const Real&t){
 Real ds=1+s,dt=1+t,ws=weighted(s);Weight wt=weight_jet(t);Shape a=shape(s,t);
 Real r2=sq(a.r),j2=sq(a.k),j3=cube(a.k)/(1-r2),j4=sq(sq(a.k))*(3-r2)/(3*cube(1-r2));
 Real numerator=ws*wt.w-sq(wt.w),mass=(numerator*j2-ws*wt.w/sq(ds))/dt;
 Real kinetic=(ws*(wt.first+t*wt.second/3)*(j2-1/sq(ds))+numerator*(-j3/3-2*j4/3)+ws*wt.first/3*((s-t+1)*j3-j2))/dt;
 return mass+kinetic;
}
Real forest(const Real&s,const Real&t,const Real&b,bool asin_kernel){
 Real kp=rt(t)*b,u=t+1-2*kp,dot=t-kp,ws=weighted(s),wt=weighted(t),wu=weighted(u);
 return (pair(t,u,dot,t*(1-sq(b)),s,asin_kernel)*(ws*wu-sq(wt))-ws*wu/sq(1+s))/(1+t);
}
Real orbit(const Real&s,const Real&t,const Real&u,bool asin_kernel){
 Real total(0);
 for(long i=0;i<8;++i){
  Real fraction=u/2+Real(i)/8;
  if(arb_equal(fraction.v,Real(1).v)){std::cerr<<"dynamic-pair-zero-weight\n";continue;}
  Real b,si;arb_cos_pi(b.v,fraction.v,precision);arb_sin_pi(si.v,fraction.v,precision);
  total=total+Real(8)/3*sq(sq(si))*(forest(s,t,b,asin_kernel)+forest(s,t,-b,asin_kernel))/2;
 }
 return total/8-subtractions(s,t);
}
Real density(const std::array<Real,3>&p,bool asin_kernel,bool diagonal){
 Real dx=1-p[0],dy=1-p[1],a=p[0]/cube(dx),b=p[1]/cube(dy),da=(1+2*p[0])/sq(sq(dx)),db=(1+2*p[1])/sq(sq(dy)),value;
 if(!diagonal)value=da*db*sq(a)*sq(b)*orbit(a,b,p[2],asin_kernel)/4;
 else{Real s=sq(a+b),t=sq(a),jac=4*a*(a+b)*da*db;value=jac*sq(s)*sq(t)*(orbit(s,t,p[2],asin_kernel)+orbit(t,s,p[2],asin_kernel))/4;}
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete original G8 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G8 recipe required");}
}
int main(int argc,char**argv){try{using namespace weighted_pair;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);std::string mode;unsigned bits,count;
 if(!(std::cin>>mode>>bits>>count)||(mode!="dynamic-pair-cubic-atan"&&mode!="dynamic-pair-cubic-asin"&&mode!="dynamic-pair-diagonal-atan"&&mode!="dynamic-pair-diagonal-asin")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded complete G8 request required");
 precision=bits;bool asin_kernel=mode.find("asin")!=std::string::npos,diagonal=mode.find("diagonal")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,3>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open three-cube required");}
 if(i)std::cout<<',';interval(density(p,asin_kernel,diagonal));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing G8 request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
