// Complete G48 forest after the full triangle parameter has been integrated.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace overlap_radial {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
#include "shared_overlap_identity.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real log(const Real&a){Real r;arb_log(r.v,a.v,precision);return r;}
Real lp(const Real&a){Real r;arb_log1p(r.v,a.v,precision);return r;}
Real cube(const Real&a){return a*sq(a);}
void pos(const Real&a,const char*msg){if(!arb_is_positive(a.v))throw std::runtime_error(msg);}
Real f0(const Real&t){pos(t,"positive root factor");return t*log(t)-t;}
Real f1(const Real&t){pos(t,"positive root factor");return sq(t)*log(t)/2-sq(t)/4;}
Real weighted_log(const Real&a,const Real&b,const Real&c,bool center){
 pos(a,"positive quadratic norm");pos(c,"positive first endpoint");Real end=c+b-a;pos(end,"positive second endpoint");
 Real disc=rt(sq(b)+4*a*c);
 if(!center){Real plus,minus;
   if(arb_is_negative(b.v)){minus=(b-disc)/(2*a);plus=-c/(a*minus);}
   else{plus=(b+disc)/(2*a);minus=-c/(a*plus);}
   pos(plus-1,"positive upper root gap");pos(-minus,"negative lower root");
   return log(a)/2+(1-plus)*(f0(plus)-f0(plus-1))+f1(plus)-f1(plus-1)
       +(1-minus)*(f0(1-minus)-f0(-minus))-f1(1-minus)+f1(-minus);
 }
 Real middle=b/(2*a),radius=disc/(2*a),w0=-middle,w1=1-middle;
 auto primitive=[&](const Real&w,const Real&p){pos(radius+w,"positive centered root");pos(radius-w,"positive centered root");
   return w*log(p)-2*w+radius*(log(radius+w)-log(radius-w));};
 // Independent centered-quadratic primitive of z*log(P), from dP=-2*a*w*dw.
 Real first=primitive(w1,end)-primitive(w0,c);
 Real second=-(end*log(end)-end-c*log(c)+c)/(2*a);
 return (1-middle)*first-second;
}
Real triangle(const Real&s,const Real&h,const Real&v,bool center){
 return -weighted_log(s+sq(v)-2*v*h,s+v*(2*v-1)-2*v*h,1+v*(1-v),center);
}
Real density(const std::array<Real,3>&p,bool center,bool square){
 Real x=square?sq(p[0]):p[0],ds=(square?2*p[0]:Real(1))/sq(1-x),s=x/(1-x),v=p[2],D=1+s;
 Real pi;arb_const_pi(pi.v,precision);Real angle=pi*p[1]/2,sine,cosine;
 arb_sin_cos(sine.v,cosine.v,angle.v,precision);Real h=rt(s)*cosine,H;arb_asinh(H.v,(rt(s)/2).v,precision);
 Real B=(s+4)*rt((s+4)/s)*H/3-4*s/9-Real(4)/3;
 Real L0=rt((s+4)/s)*H-1,J0=4*H/rt(s*(s+4)),J00=2*(1+J0)/(s+4);
 Real I20=J0/2-(1-J0)/s,I30=J0/2-3*(1-J0)/(2*s);
 Real L2=v*I20-sq(v)*I30-sq(h)*sq(v)/sq(s)*(1-2*J0+J00);
 Real T0=-L0/D,T2=-L2/D+2*sq(h)*v*(1-J0)/(s*sq(D))-L0*(4*sq(h)/cube(D)-1/sq(D));
 Real paired=(triangle(s,h,v,center)/(D+1-2*h)+triangle(s,-h,v,center)/(D+1+2*h))/2;
 Real sigma=v*(1-v),counter=((1+sigma)*lp(sigma)-sigma)/sq(D);
 Real value=ds*sq(s)*sq(B)/(2*cube(D))*Real(8)/3*sq(sq(sine))*(paired-T0-T2-counter);
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete overlap remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G48 recipe required");}
}
int main(int argc,char**argv){try{using namespace overlap_radial;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||
 (mode!="overlap-linear-root"&&mode!="overlap-linear-center"&&mode!="overlap-square-root"&&mode!="overlap-square-center")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded three-dimensional overlap request required");
 precision=bits;bool center=mode.find("center")!=std::string::npos,square=mode.find("square")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,3>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open three-cube required");}
 if(i)std::cout<<',';interval(density(p,center,square));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing overlap request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
