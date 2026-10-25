// Complete G48 forest after the full triangle parameter has been integrated.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace mixed_double_triangle {
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
Real control(const Real&s,const Real&c,const Real&v){
 Real sigma=v*(1-v),ell=lp(sigma),c2=sq(c),c4=sq(c2),v2=sq(v),v3=v*v2,v4=sq(v2);
 Real a=-(1-12*c2+16*c4)/2-sigma*((v2+1)*(4*c2-1)+4*v*c2)-sq(sigma)/2;
 Real b=c4*(2*v4+Real(8)*v3/3+4*v2+8*v+16)+c2*((-4*v4+4*v2+8*v+4)*ell-8*v4+v2-4*v-12)
       +(Real(3)*v4/2-2*v3-v2/2-Real(1)/2)*ell+v4/2+v3/2+v2/4-v/2+1;
 Real half_asinh;arb_asinh(half_asinh.v,(Real(1)/2).v,precision);Real mean=-Real(16)/27+5*rt(Real(5))*half_asinh/9;
 return sq(s)/(cube(1+s)*sq(1+s))*(a*lp(s)+b-mean);
}
Real one(const Real&s,const Real&c,const Real&v,const Real&w,bool center){
 Real h=rt(s)*c,d=1+s,z;arb_asinh(z.v,(rt(s)/2).v,precision);
 Real l0=rt((s+4)/s)*z-1,j=4*z/rt(s*(s+4)),j2=2*(1+j)/(s+4),a=(1-j)/s;
 Real i20=j/2-a,i30=j/2-3*a/2;
 auto l2=[&](const Real&u){return u*i20-sq(u)*i30-sq(h)*sq(u)/sq(s)*(1-2*j+j2);};
 Real c0=sq(l0),c1=-l0*h*(v+w)*a,c2=l0*(l2(v)+l2(w))+sq(h)*v*w*sq(a);
 Real mass=c0/d,kinetic=c2/d+2*h*c1/sq(d)+c0*(4*sq(h)/cube(d)-1/sq(d));
 Real paired=(triangle(s,h,v,center)*triangle(s,h,w,center)/(d+1-2*h)+triangle(s,-h,v,center)*triangle(s,-h,w,center)/(d+1+2*h))/2;
 auto bp=[&](const Real&u){Real sigma=u*(1-u);return (1+sigma)*lp(sigma)-sigma;};
 Real parent=-l0/sq(d)*(bp(v)+bp(w)),zero_mean=-l0*(control(s,c,v)+control(s,c,w));
 return paired-mass-kinetic-parent-zero_mean;
}
Real density(const std::array<Real,4>&p,bool center,bool eight){
 Real d=1-p[0],s=p[0]/cube(d),ds=(1+2*p[0])/sq(sq(d)),total(0);unsigned count=eight?8:4;
 for(unsigned i=0;i<count;++i){
  Real fraction=p[1]/2+Real(long(i))/Real(long(count));
  if(arb_equal(fraction.v,Real(1).v)){std::cerr<<"mixed-double-triangle-zero-weight\n";continue;}
  Real c,si;arb_cos_pi(c.v,fraction.v,precision);arb_sin_pi(si.v,fraction.v,precision);
  Real local=one(s,c,p[2],p[3],center)+one(s,c,1-p[2],p[3],center)+one(s,c,p[2],1-p[3],center)+one(s,c,1-p[2],1-p[3],center);
  total=total+Real(8)/3*sq(sq(si))*local/4;
 }
 Real value=ds*sq(s)*bubble(s)/(2*sq(1+s))*total/Real(long(count));
 if(!arb_is_finite(value.v))throw std::runtime_error("finite complete original G63 remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G63 recipe required");}
}
int main(int argc,char**argv){try{using namespace mixed_double_triangle;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||
 (mode!="mixed-double-triangle-four-root"&&mode!="mixed-double-triangle-four-center"&&mode!="mixed-double-triangle-eight-root"&&mode!="mixed-double-triangle-eight-center")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded four-dimensional double-triangle request required");
 precision=bits;bool center=mode.find("center")!=std::string::npos,eight=mode.find("eight")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,4>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open four-cube required");}
 if(i)std::cout<<',';interval(density(p,center,eight));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing box request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
