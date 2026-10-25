// Complete G48 forest after the full triangle parameter has been integrated.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace dynamic_triangle {
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


Real bubble(const Real&s){Real h;arb_asinh(h.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*h/3-4*s/9-Real(4)/3;}
Real weighted(const Real&s){return bubble(s)/sq(1+s);}
Real one(const Real&s,const Real&h,const Real&v,bool center){
 Real d=1+s,asinh;arb_asinh(asinh.v,(rt(s)/2).v,precision);Real b=bubble(s),j0=4*asinh/rt(s*(s+4)),j00=2*(1+j0)/(s+4);
 Real bp=rt(s+4)*(s-2)*asinh/(3*s*rt(s))+(s+4)/(6*s)-Real(4)/9,bpp=1/(6*s)-(1-j0)/sq(s);
 Real w=b/sq(d),wp=bp/sq(d)-2*b/cube(d),wpp=bpp/sq(d)-4*bp/cube(d)+6*b/sq(sq(d));
 Real l0=rt((s+4)/s)*asinh-1,i20=j0/2-(1-j0)/s,i30=j0/2-3*(1-j0)/(2*s);
 Real l2=v*i20-sq(v)*i30-sq(h)*sq(v)/sq(s)*(1-2*j0+j00);
 Real mass=-w*l0,kinetic=-w*l2-2*sq(h)*v*(1-j0)/s*wp-l0*(wp+2*sq(h)*wpp);
 Real paired=(weighted(s+1-2*h)*triangle(s,h,v,center)+weighted(s+1+2*h)*triangle(s,-h,v,center))/2;
 Real sigma=v*(1-v),counter=w/d*((1+sigma)*lp(sigma)-sigma);
 return paired-mass-kinetic-counter;
}
Real density(const std::array<Real,3>&p,bool center,bool eight){
 Real d=1-p[0],s=p[0]/cube(d),ds=(1+2*p[0])/sq(sq(d)),total(0);unsigned count=eight?8:4;
 for(unsigned i=0;i<count;++i){
  Real fraction=p[1]/2+Real(long(i))/Real(long(count));
  if(arb_equal(fraction.v,Real(1).v)){std::cerr<<"dynamic-triangle-zero-weight\n";continue;}
  Real co,si;arb_cos_pi(co.v,fraction.v,precision);arb_sin_pi(si.v,fraction.v,precision);Real h=rt(s)*co;
  total=total+Real(8)/3*sq(sq(si))*(one(s,h,p[2],center)+one(s,h,1-p[2],center))/2;
 }
 Real value=ds*sq(s)*weighted(s)*total/(2*Real(long(count)));
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete shifted-weight G17 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G17 recipe required");}
}
int main(int argc,char**argv){try{using namespace dynamic_triangle;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||
 (mode!="dynamic-triangle-four-root"&&mode!="dynamic-triangle-four-center"&&mode!="dynamic-triangle-eight-root"&&mode!="dynamic-triangle-eight-center")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded three-dimensional box request required");
 precision=bits;bool center=mode.find("center")!=std::string::npos,eight=mode.find("eight")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,3>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open three-cube required");}
 if(i)std::cout<<',';interval(density(p,center,eight));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing box request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
