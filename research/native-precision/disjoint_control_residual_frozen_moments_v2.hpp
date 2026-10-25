// Analytically integrate the last simplex parameter, including a rigorous tail.
// Reuse the frozen graph recipe and bubble implementation without editing it.
#define main disjoint_five_dimensional_main
#include "disjoint-triangle-bubble-native-v4/integrate.cpp"
#undef main
namespace disjoint_v_integrated {
using arb_cw::Real;using arb_cw::precision;
using disjoint_interval::Recipe;using fourier_e2_root::Z;
using fourier_e2_root::readq;using fourier_e2_root::interval;
Real square(const Real&x){Real r;arb_sqr(r.v,x.v,precision);return r;}
Real power(const Real&x,unsigned n){Real r;arb_pow_ui(r.v,x.v,n,precision);return r;}
Real exponential(const Real&x){Real r;arb_exp(r.v,x.v,precision);return r;}
Real expm1(const Real&x){Real r;arb_expm1(r.v,x.v,precision);return r;}
Real log_one_plus(const Real&x){Real r;arb_log1p(r.v,x.v,precision);return r;}
Real logarithm(const Real&x){Real r;arb_log(r.v,x.v,precision);return r;}
Real root(const Real&x){Real r;arb_sqrt(r.v,x.v,precision);return r;}
Real absolute(const Real&x){Real r;arb_abs(r.v,x.v);return r;}
void positive(const Real&x,const char*label){if(!arb_is_positive(x.v))throw std::runtime_error(label);}
struct Moment{Real beta,gamma,rho,tail,sum;};
Moment moments(const Real&fc,const Real&slope,const Real&curvature){
 positive(fc,"positive centered triangle mass");Moment out;
 out.beta=slope/(2*fc);out.gamma=curvature/(4*fc);
 if(!arb_is_nonnegative(out.gamma.v))throw std::runtime_error("nonnegative quadratic curvature");
 out.rho=(root(square(out.beta)+4*out.gamma)+absolute(out.beta))/2;
 if(!arb_lt(out.rho.v,(Real(2)/3).v))throw std::runtime_error("unit-leg quadratic moment bound");
 Real coefficient=square(out.beta)+2*out.gamma,prior(2),current=coefficient;
 for(unsigned n=1;n<=128;++n){
  out.sum=out.sum+current/Real(long(2*n*(2*n+1)));
  Real next=coefficient*current-square(out.gamma)*prior;prior=current;current=next;
 }
 out.tail=2*power(out.rho,258)/(Real(258*259)*(1-square(out.rho)));
 Real remainder;arb_union(remainder.v,Real(0).v,out.tail.v,precision);out.sum=out.sum+remainder;
 return out;
}
struct Value{Real value;std::array<Moment,2>moments;};
Value density(const Recipe&r,unsigned chart,const std::array<Real,4>&p){
 Real t=exponential(p[0]*r.ln2),z=p[1],u=p[3],y,jac;
 positive(t,"positive t");
 if(chart==0){Real L=log_one_plus(t/2);y=expm1(u*L)/t;jac=L*exponential(u*L)/t;}
 else{y=(1+u)/2;jac=Real(1)/2;}
 Real co,si;arb_sin_cos_pi(si.v,co.v,p[2].v,precision);Real second=(-co+r.sqrt3*si)/2;
 std::array<Real,2>projection;std::array<long,2>norms;long inside=0;
 for(const auto&a:r.pairs){
  if(a.i==0){projection[a.j-1]=a.c*(a.a*co+a.b*second);norms[a.j-1]=a.norm;}
  else inside=a.norm;
 }
 if(inside!=1)throw std::invalid_argument("unit triangle inner external leg");
 std::array<std::array<long,2>,2>shift;
 for(const auto&a:r.pairs)if(a.i==0)shift[a.j-1]={a.c*a.a,a.c*a.b};
 long sa=shift[0][0]-shift[1][0],sb=shift[0][1]-shift[1][1];
 if(sa*sa-sa*sb+sb*sb!=1)throw std::invalid_argument("coherent unit triangle shift difference");
 Real alpha=y*(1-y),curvature=square(1-y),M=1+t*alpha;
 Real Cc=alpha*Real(norms[0]+norms[1])/2+curvature/4;
 Real Cd=alpha*Real(norms[0]-norms[1]);
 Real q=root(t*z),Lc=q*alpha*(projection[0]+projection[1]);
 Real Ld=2*q*alpha*(projection[0]-projection[1]);
 Value result;result.moments[0]=moments(M+Cc+Lc,Cd+Ld,curvature);
 result.moments[1]=moments(M+Cc-Lc,Cd-Ld,curvature);
 Real c=Cc/M,l=Lc/M;
 Real even=log_one_plus(2*c+square(c)-square(l))/2-(result.moments[0].sum+result.moments[1].sum)/2;
 Real arg=Lc/(M+Cc),odd;arb_atanh(odd.v,arg.v,precision);
 odd=odd-(result.moments[0].sum-result.moments[1].sum)/2;
 Real h=r.host[0][0]*co+r.host[0][1]*second,k=r.host[1][0]*co+r.host[1][1]*second;
 Real A=1+t,X=2+t,invX=1/X,eta=4*(t*invX)*invX*z,beta=(2*A+1)*square(invX);
 Real da=1-eta*square(h),db=1-eta*square(k);positive(da,"positive normalized host A");positive(db,"positive normalized host B");
 Real U=1+eta*h*k,Dnum=eta*square(h+k)-U*(beta+eta*h*k),AX=A*invX;
 Z bc=disjoint_interval::bubble(Z(t));Real B;arb_set(B.v,acb_realref(bc.v));
 if(!arb_contains_zero(acb_imagref(bc.v)))throw std::runtime_error("real bubble check");
 Real prefactor=power(t/A,3)*(B/A)/da/db;
 Real odd_coefficient=-square(AX)*(2*q*invX)*(h+k);
 result.value=-r.ln2*(1-z)*(1-y)*jac*prefactor*(Dnum*logarithm(M)+square(AX)*U*even+odd_coefficient*odd);
 if(!arb_is_finite(result.value.v))throw std::runtime_error("finite analytic v integral");
 return result;
}
int execute(int argc,char**argv){
 if(argc!=2)throw std::invalid_argument("recipe path required");std::string mode;unsigned bits,chart;
 if(!(std::cin>>mode>>bits>>chart)||mode!="point4"||bits<256||bits>2048||chart>1)throw std::invalid_argument("point4 precision and chart");
 precision=bits;Recipe r;disjoint_interval::prepare(r,argv[1]);std::array<Real,4>p;
 for(unsigned j=0;j<4;++j){p[j]=readq();if(!arb_ge(p[j].v,Real(j?0:-12).v)||!arb_le(p[j].v,Real(j?1:64).v))throw std::invalid_argument("physical four-dimensional point");}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing point4 request");
 Value out=density(r,chart,p);
 std::cout<<"{\"accepted\":true,\"mode\":\"point4\",\"bits\":"<<bits<<",\"chart\":"<<chart<<",\"calls\":1,\"moment_order\":128,\"value\":";
 disjoint_interval::print_z(Z(out.value));std::cout<<",\"moments\":[";
 for(unsigned j=0;j<2;++j){if(j)std::cout<<',';const auto&m=out.moments[j];
  std::cout<<"{\"beta\":";interval(m.beta);std::cout<<",\"gamma\":";interval(m.gamma);std::cout<<",\"rho\":";interval(m.rho);std::cout<<",\"tail\":";interval(m.tail);std::cout<<",\"sum\":";interval(m.sum);std::cout<<'}';
 }
 std::cout<<"]}\n";return 0;
}
}
