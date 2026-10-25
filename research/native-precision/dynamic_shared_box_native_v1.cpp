// Complete G48 forest after the full triangle parameter has been integrated.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace dynamic_box {
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

Real box(const Real&s,const Real&h,const Real&v,bool center){
 if(!arb_contains_zero(h.v))return (triangle(s,h,v,center)-triangle(s,-h,v,center))/h;
 Real a=s+sq(v),b=s+v*(2*v-1),c=1+v*(1-v),delta=rt(sq(b)+4*a*c),x=(2*a-b)/delta,y=b/delta,q0;
 if(center)q0=(real_log((1+x)/(1-x))+real_log((1+y)/(1-y)))/delta;
 else {Real left,right;arb_atanh(left.v,x.v,precision);arb_atanh(right.v,y.v,precision);q0=2*(left+right)/delta;}
 Real q1=(b*q0-real_log((c+b-a)/c))/(2*a),q2=(c*q0+b*q1-1)/a,q3=(c*q1+b*q2-Real(1)/2)/a;
 Real bound=2*cube(v)*sq(h)/105,tail;arb_union(tail.v,(-bound).v,bound.v,precision);
 std::cerr<<"dynamic-box-removable-limit\n";return 4*v*(q1-2*q2+q3)+tail;
}
Real one(const Real&s,const Real&h,const Real&v,bool center){
 Real root=rt(s),asinh;arb_asinh(asinh.v,(root/2).v,precision);
 Real j1=4*asinh/rt(s*(s+4)),j2=2*(1+j1)/(s+4),j3=(1+3*j2)/(s+4);
 Real m=(1-j1)/s,i1=(j1-j2)/s,i2=(1-2*j1+j2)/sq(s),i3=(1-3*j1+3*j2-j3)/cube(s);
 Real mass=2*v*m,kinetic=-2*sq(v)*(i1-2*i2)+2*cube(v)*(i1-3*i2)+Real(8)/3*sq(h)*cube(v)*i3;
 Real sigma=v*(1-v),counter=2*((1+sigma)*lp(sigma)-sigma)/sq(1+s);
 return box(s,h,v,center)-mass-kinetic-counter;
}
Real density(const std::array<Real,3>&p,bool center,bool eight){
 Real d=1-p[0],s=p[0]/cube(d),ds=(1+2*p[0])/sq(sq(d)),total(0);unsigned count=eight?8:4;
 for(unsigned i=0;i<count;++i){
  Real fraction=p[1]/2+Real(long(i))/Real(long(count));
  if(arb_equal(fraction.v,Real(1).v)){std::cerr<<"dynamic-box-zero-weight\n";continue;}
  Real co,si;arb_cos_pi(co.v,fraction.v,precision);arb_sin_pi(si.v,fraction.v,precision);Real h=rt(s)*co;
  total=total+Real(8)/3*sq(sq(si))*(one(s,h,p[2],center)+one(s,h,1-p[2],center))/2;
 }
 Real asinh;arb_asinh(asinh.v,(rt(s)/2).v,precision);Real B=(s+4)*rt((s+4)/s)*asinh/3-4*s/9-Real(4)/3;
 Real value=ds*sq(s)*sq(B)/(2*cube(1+s))*total/Real(long(count));
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete G1 box remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G1 recipe required");}
}
int main(int argc,char**argv){try{using namespace dynamic_box;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||
 (mode!="dynamic-box-four-root"&&mode!="dynamic-box-four-center"&&mode!="dynamic-box-eight-root"&&mode!="dynamic-box-eight-center")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded three-dimensional box request required");
 precision=bits;bool center=mode.find("center")!=std::string::npos,eight=mode.find("eight")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,3>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open three-cube required");}
 if(i)std::cout<<',';interval(density(p,center,eight));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing box request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
