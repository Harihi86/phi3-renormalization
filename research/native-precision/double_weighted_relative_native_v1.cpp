// G20 complete original Taylor remainder after exact transverse angular integration.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace relative_e2 {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
#include "shared_overlap_identity.hpp"
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real absolute(const Real&a){Real r;arb_abs(r.v,a.v);return r;}
Real cube(const Real&a){return a*sq(a);}
std::array<Real,3> transverse(const Real&c,const Real&w,bool atanh_mode){
 if(!arb_is_positive((c-w).v))throw std::runtime_error("positive massive transverse denominator required");
 if(arb_is_zero(w.v))return {Real(1)/c,-Real(1)/sq(c),Real(2)/cube(c)};
 Real h;
 if(atanh_mode)arb_atanh(h.v,(w/c).v,precision);
 else{Real plus,minus;arb_log1p(plus.v,(w/c).v,precision);arb_log1p(minus.v,(-w/c).v,precision);h=(plus-minus)/2;}
 Real value=Real(3)/(2*sq(w))*(c-(sq(c)-sq(w))*h/w);
 Real first=Real(3)/sq(w)*(1-c*h/w),second=Real(3)/sq(w)*(c/(sq(c)-sq(w))-h/w);
 return {value,first,second};
}
Real remainder(const Real&s,const Real&t,const Real&a,const Real&b,const Real&width,bool mode){
 Real ds=1+s,dt=1+t,rs=rt(s),rr=rt(t),la=-2*rs*a,lb=2*rr*b;
 Real c=1+s+t-2*rs*rr*a*b,linear=la+lb;auto f=transverse(c,width,mode);
 Real g0=Real(1)/(ds*dt),g1=-(la/ds+lb/dt)/(ds*dt);
 Real g2=(sq(la)/sq(ds)-1/ds+sq(lb)/sq(dt)-1/dt+la*lb/(ds*dt))/(ds*dt);
 Real mass=g0*f[0],kinetic=g2*f[0]+g1*linear*f[1]+g0*(f[1]+sq(linear)*f[2]/2);
 Real full=(transverse(c+linear+1,width,mode)[0]/((ds+la+1)*(dt+lb+1))
           +transverse(c-linear+1,width,mode)[0]/((ds-la+1)*(dt-lb+1)))/2;
 Real ps=((Real(1)/(ds+la+1)+Real(1)/(ds-la+1))/2-1/ds-(sq(la)/cube(ds)-1/sq(ds)))/sq(dt);
 Real pt=((Real(1)/(dt+lb+1)+Real(1)/(dt-lb+1))/2-1/dt-(sq(lb)/cube(dt)-1/sq(dt)))/sq(ds);
 return full-mass-kinetic-ps-pt;
}
Real bubble(const Real&s){Real h;arb_asinh(h.v,(rt(s)/2).v,precision);return (s+4)*rt((s+4)/s)*h/3-4*s/9-Real(4)/3;}
Real density(const std::array<Real,4>&p,bool mode,bool eight){
 Real dx=1-p[0],dy=1-p[1],s=p[0]/cube(dx),t=p[1]/cube(dy);
 Real ds=(1+2*p[0])/sq(sq(dx)),dt=(1+2*p[1])/sq(sq(dy));
 Real pi;arb_const_pi(pi.v,precision);unsigned count=eight?8:4;Real total(0);
 for(unsigned ia=0;ia<count;++ia){
  Real fa=p[2]+Real(long(2*ia))/Real(long(count));if(arb_equal(fa.v,Real(1).v)||arb_equal(fa.v,Real(2).v))continue;
  Real sa,a;arb_sin_cos(sa.v,a.v,(pi*fa).v,precision);Real wa=Real(8)/3*sq(sq(sa));
  for(unsigned ib=0;ib<count;++ib){
   Real fb=p[3]+Real(long(2*ib))/Real(long(count));if(arb_equal(fb.v,Real(1).v)||arb_equal(fb.v,Real(2).v))continue;
   Real sb,b;arb_sin_cos(sb.v,b.v,(pi*fb).v,precision);Real wb=Real(8)/3*sq(sq(sb));
   Real width=2*rt(s*t)*absolute(sa*sb);total=total+wa*wb*remainder(s,t,a,b,width,mode);
  }
 }
 Real value=ds*dt*sq(s)*sq(t)*bubble(s)*bubble(t)/(4*sq(1+s)*sq(1+t))*total/Real(long(count*count));
 if(!arb_is_finite(value.v))throw std::runtime_error("finite full relative-angle remainder required");return value;
}
void recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra)||id!=GRAPH_ID||sha!=GRAPH_SHA)throw std::invalid_argument("exact original G20 recipe required");}
}
int main(int argc,char**argv){try{using namespace relative_e2;
 if(argc!=2)throw std::invalid_argument("one recipe required");recipe(argv[1]);
 std::string mode;unsigned bits,count;
 if(!(std::cin>>mode>>bits>>count)||
 (mode!="relative-four-log"&&mode!="relative-four-atanh"&&mode!="relative-eight-log"&&mode!="relative-eight-atanh")||bits<512||bits>2048||count<1||count>512)
 throw std::invalid_argument("bounded four-dimensional relative-angle request required");
 precision=bits;bool atanh_mode=mode.find("atanh")!=std::string::npos,eight=mode.find("eight")!=std::string::npos;
 std::cout<<"{\"id\":\""<<GRAPH_ID<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,4>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open four-cube required");}
 if(i)std::cout<<',';interval(density(p,atanh_mode,eight));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing relative-angle request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
