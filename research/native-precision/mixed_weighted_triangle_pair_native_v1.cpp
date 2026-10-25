// Original G46: full spherical pair, complete triangle, both parent forests and T2.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace mixed_triangle_pair {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
#include "shared_overlap_identity.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real real_log(const Real&a){Real r;arb_log(r.v,a.v,precision);return r;}
Real lp(const Real&a){Real r;arb_log1p(r.v,a.v,precision);return r;}
Real cube(const Real&a){return a*sq(a);}
void pos(const Real&a,const char*msg){if(!arb_is_positive(a.v))throw std::runtime_error(msg);}
Real f0(const Real&t){pos(t,"positive root factor");return t*real_log(t)-t;}
Real f1(const Real&t){pos(t,"positive root factor");return sq(t)*real_log(t)/2-sq(t)/4;}
Real weighted_log(const Real&a,const Real&b,const Real&c,bool center){
 pos(a,"positive quadratic norm");pos(c,"positive first endpoint");Real end=c+b-a;pos(end,"positive second endpoint");
 Real disc=rt(sq(b)+4*a*c);
 if(!center){Real plus,minus;
   if(arb_is_negative(b.v)){minus=(b-disc)/(2*a);plus=-c/(a*minus);}
   else{plus=(b+disc)/(2*a);minus=-c/(a*plus);}
   pos(plus-1,"positive upper root gap");pos(-minus,"negative lower root");
   return real_log(a)/2+(1-plus)*(f0(plus)-f0(plus-1))+f1(plus)-f1(plus-1)
       +(1-minus)*(f0(1-minus)-f0(-minus))-f1(1-minus)+f1(-minus);
 }
 Real middle=b/(2*a),radius=disc/(2*a),w0=-middle,w1=1-middle;
 auto primitive=[&](const Real&w,const Real&p){pos(radius+w,"positive centered root");pos(radius-w,"positive centered root");
   return w*real_log(p)-2*w+radius*(real_log(radius+w)-real_log(radius-w));};
 // Independent centered-quadratic primitive of z*real_log(P), from dP=-2*a*w*dw.
 Real first=primitive(w1,end)-primitive(w0,c);
 Real second=-(end*real_log(end)-end-c*real_log(c)+c)/(2*a);
 return (1-middle)*first-second;
}
Real triangle(const Real&s,const Real&h,const Real&v,bool center){
 return -weighted_log(s+sq(v)-2*v*h,s+v*(2*v-1)-2*v*h,1+v*(1-v),center);
}



Real ab(const Real&a){Real r;arb_abs(r.v,a.v);return r;}
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

Real bubble(const Real&s){Real z;arb_asinh(z.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*z/3-4*s/9-Real(4)/3;}
Real weighted(const Real&s){return bubble(s)/sq(1+s);}
std::array<Real,3> triangle_jet(const Real&s,const Real&h,const Real&v){
 Real z;arb_asinh(z.v,(rt(s)/2).v,precision);Real l0=rt((s+4)/s)*z-1,j=4*z/rt(s*(s+4)),j2=2*(1+j)/(s+4),a=(1-j)/s;
 Real i20=j/2-a,i30=j/2-3*a/2,l2=v*i20-sq(v)*i30-sq(h)*sq(v)/sq(s)*(1-2*j+j2);
 return {-l0,h*v*a,-l2};
}
std::array<Real,3> pair_jet(const Real&s,const Real&t,const Real&h){
 Shape shape0=shape(s,t);Real k=shape0.k,delta=rt(sq(1+s+t)-4*s*t),u=1+(1+s-t)/delta;
 Real first=-cube(k)*u,second=Real(3)/2*sq(sq(k))*sq(u)-4*t*cube(k)/cube(delta),j3=sq(k)/delta;
 return {sq(k)-1/sq(1+t),-h*first,(first+2*sq(h)*second-j3)/3};
}
Real bubble_radial(const Real&s){
 Real d=1+s;Shape a=shape(s,Real(1));Real mean=a.k*(1-sq(a.r)/3),kin=-1/sq(d)+2*s/(3*cube(d));
 return (mean-1/d-kin)/d;
}
Real one(const Real&s,const Real&t,const Real&c,const Real&v,bool center){
 Real h=rt(s)*c,d=1+s;auto a=pair_jet(s,t,h),b=triangle_jet(s,h,v);
 Real c0=a[0]*b[0],c1=a[0]*b[1]+a[1]*b[0],c2=a[0]*b[2]+a[1]*b[1]+a[2]*b[0];
 Real mass=c0/sq(d),kin=c2/sq(d)+2*h*c1/cube(d)+c0*(4*sq(h)/sq(sq(d))-1/cube(d));
 auto product=[&](const Real&sign){
  Real hh=sign*h,u=s+1-2*hh,dot=s-hh;
  Real finite=pair(s,u,dot,s*(1-sq(c)),t,center)-1/sq(1+t);
  return finite*triangle(s,hh,v,center)/(d*(d+1-2*hh));
 };
 Real sigma=v*(1-v),bp=(1+sigma)*lp(sigma)-sigma;
 Real parent3=a[0]/cube(d)*bp,parent4=triangle_jet(t,Real(0),v)[0]/sq(1+t)*bubble_radial(s);
 return (product(Real(1))+product(Real(-1)))/2-mass-kin-parent3-parent4;
}
Real orbit(const Real&s,const Real&t,const Real&u,const Real&v,bool center){
 Real total(0);
 for(long i=0;i<8;++i){
  Real fraction=u/2+Real(i)/8;
  if(arb_equal(fraction.v,Real(1).v)){std::cerr<<"mixed-triangle-pair-zero-weight\n";continue;}
  Real c,si;arb_cos_pi(c.v,fraction.v,precision);arb_sin_pi(si.v,fraction.v,precision);
  total=total+Real(8)/3*sq(sq(si))*(one(s,t,c,v,center)+one(s,t,c,1-v,center))/2;
 }
 return total/8;
}
Real density(const std::array<Real,4>&p,bool center,bool diagonal){
 Real dx=1-p[0],dy=1-p[1],a=p[0]/cube(dx),b=p[1]/cube(dy),da=(1+2*p[0])/sq(sq(dx)),db=(1+2*p[1])/sq(sq(dy)),value;
 if(!diagonal)value=da*db*sq(a)*sq(b)*weighted(b)*orbit(a,b,p[2],p[3],center)/4;
 else{Real s=sq(a+b),t=sq(a),jac=4*a*(a+b)*da*db;value=jac*sq(s)*sq(t)*(weighted(t)*orbit(s,t,p[2],p[3],center)+weighted(s)*orbit(t,s,p[2],p[3],center))/4;}
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete original G46 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G46 recipe required");}
}
int main(int argc,char**argv){try{using namespace mixed_triangle_pair;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||
 (mode!="mixed-triangle-pair-cubic-root"&&mode!="mixed-triangle-pair-cubic-center"&&mode!="mixed-triangle-pair-diagonal-root"&&mode!="mixed-triangle-pair-diagonal-center")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded four-dimensional pair-triangle request required");
 precision=bits;bool center=mode.find("center")!=std::string::npos,diagonal=mode.find("diagonal")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,4>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open four-cube required");}
 if(i)std::cout<<',';interval(density(p,center,diagonal));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing box request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
