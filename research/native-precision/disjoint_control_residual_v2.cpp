// Exact frozen mathematical sections; separate batch entry point.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace disjoint_original_central {
using arb_cw::Real;using arb_cw::precision;
using disjoint_interval::Recipe;using disjoint_interval::Rule;
using fourier_e2_root::readq;using fourier_e2_root::interval;
Real sq(const Real&x){Real r;arb_sqr(r.v,x.v,precision);return r;}
Real pw(const Real&x,unsigned n){Real r;arb_pow_ui(r.v,x.v,n,precision);return r;}
Real ex(const Real&x){Real r;arb_exp(r.v,x.v,precision);return r;}
Real em1(const Real&x){Real r;arb_expm1(r.v,x.v,precision);return r;}
Real lp(const Real&x){Real r;arb_log1p(r.v,x.v,precision);return r;}
Real lg(const Real&x){Real r;arb_log(r.v,x.v,precision);return r;}
Real rt(const Real&x){Real r;arb_sqrt(r.v,x.v,precision);return r;}
Real ah(const Real&x){Real r;arb_atanh(r.v,x.v,precision);return r;}
void pos(const Real&x,const char*what){if(!arb_is_positive(x.v))throw std::runtime_error(what);}
void validate(const Recipe&r){
 std::array<std::array<long,2>,2>shifts{};long inside=0;
 for(const auto&p:r.pairs){if(p.i==0)shifts[p.j-1]={p.c*p.a,p.c*p.b};else inside=p.norm;}
 long a=shifts[0][0]-shifts[1][0],b=shifts[0][1]-shifts[1][1];
 if(inside!=1||a*a-a*b+b*b!=1)throw std::invalid_argument("original triangle unit separation");
}
Real original_variance(const Recipe&r,const Real&t,const Real&y,const Real&q1,const Real&q2,const Real&v,long sign){
 std::array<Real,3>xs={y,(1-y)*v,(1-y)*(1-v)};Real out(1);
 for(const auto&p:r.pairs)out=out+xs[p.i]*xs[p.j]*(p.c*p.c*t+p.norm+2*sign*p.c*(p.a*q1+p.b*q2));
 pos(out,"positive original triangle mass");return out;
}
Real integrated_log(const Real&f0,const Real&fc,const Real&f1){
 // f(v)=f0+b*v-a*v^2; coefficients inferred from the original three values.
 Real a=4*fc-2*(f0+f1),b=f1-f0+a;
 if(arb_is_zero(a.v)){
  if(arb_is_zero(b.v))return lg(f0);
  if(arb_contains_zero(b.v))throw std::runtime_error("uncertain linear primitive");
  return (f1*lg(f1)-f0*lg(f0))/b-1;
 }
 pos(a,"positive inferred curvature");Real D=rt(sq(b)+4*a*f0),den=2*f0+b;
 pos(den,"positive primitive denominator");Real w=D/den;
 pos(w,"positive primitive atanh argument");pos(1-w,"primitive atanh argument below one");
 return b/(2*a)*lg(f0)+(1-b/(2*a))*lg(f1)-2+D/a*ah(w);
}
Real original_host(const Recipe&r,const Real&t,const Real&q1,const Real&q2,long sign){
 Real out=1/sq(1+t);
 for(const auto&h:r.host){long a=h[0],b=h[1];Real denominator=1+t+a*a-a*b+b*b+2*sign*(a*q1+b*q2);pos(denominator,"positive original host denominator");out=out/denominator;}
 return out;
}
Real original_bubble(const Real&t){
 Real u=rt(t/(t+4));pos(u,"positive independent bubble argument");pos(1-u,"independent bubble argument below one");
 return (t+4)/3*(ah(u)/u-1)-t/9;
}
Real density(const Recipe&r,unsigned chart,const std::array<Real,4>&p){
 Real t=ex(p[0]*r.ln2),z=p[1],u=p[3],y,jac;pos(t,"positive radial variable");
 if(chart==0){Real L=lp(t/2);y=em1(u*L)/t;jac=L*ex(u*L)/t;}else{y=(1+u)/2;jac=Real(1)/2;}
 // The endpoint is exactly zero before a singular-looking primitive is used.
 if(arb_is_zero((1-y).v)||arb_is_zero((1-z).v))return Real(0);
 Real co,si;arb_sin_cos_pi(si.v,co.v,p[2].v,precision);
 Real q=rt(t*z),q1=q*co,q2=q*(-co+r.sqrt3*si)/2,paired;
 for(long sign:{1L,-1L}){
  Real f0=original_variance(r,t,y,q1,q2,Real(0),sign);
  Real fc=original_variance(r,t,y,q1,q2,Real(1)/2,sign);
  Real f1=original_variance(r,t,y,q1,q2,Real(1),sign);
  paired=paired+original_host(r,t,q1,q2,sign)*integrated_log(f0,fc,f1)/2;
 }
 Real kernel=-original_bubble(t)*(paired-lp(y*(1-y)*t)/pw(1+t,4));
 Real value=pw(t,3)*r.ln2*(1-z)*(1-y)*jac*kernel;
 if(!arb_is_finite(value.v))throw std::runtime_error("nonfinite independent central value");return value;
}

}
namespace disjoint_control_residual {
namespace M=disjoint_v_integrated;namespace O=disjoint_original_central;
using arb_cw::Real;using arb_cw::precision;using disjoint_interval::Recipe;
using fourier_e2_root::readq;using fourier_e2_root::interval;
Real value(const Recipe&r,const std::string&mode,unsigned chart,const std::array<Real,4>&cube){
    const Real&s=cube[0],&z=cube[1],&theta=cube[2],&u=cube[3];
    Real t=M::square(s/(1-s)),x=2*(M::logarithm(s)-M::log_one_plus(-s))/r.ln2;
    std::array<Real,4>p{x,z,theta,u};Real y,jac;
    if(chart==0){Real l=M::log_one_plus(t/2);y=M::expm1(u*l)/t;jac=l*M::exponential(u*l)/t;}
    else{y=(1+u)/2;jac=Real(1)/2;}
    Real co,si;arb_sin_cos_pi(si.v,co.v,theta.v,precision);
    Real q=M::root(t*z),q1=q*co,q2=q*(-co+r.sqrt3*si)/2;
    Real A=1+t,alpha=y*(1-y),mass=1+t*alpha,point,B,full;
    if(mode=="stable"){
        Real h=r.host[0][0]*co+r.host[0][1]*(-co+r.sqrt3*si)/2;
        Real k=r.host[1][0]*co+r.host[1][1]*(-co+r.sqrt3*si)/2;
        std::array<Real,2>proj;long nsum=0;
        for(const auto&a:r.pairs)if(a.i==0){proj[a.j-1]=a.c*(a.a*co+a.b*(-co+r.sqrt3*si)/2);nsum+=a.norm;}
        Real cbar=alpha*Real(nsum)/2+M::square(1-y)/6;
        point=(4*t*z*(M::square(h)+M::square(k)+h*k)/M::square(A)-2/A)*M::logarithm(mass)+cbar/mass;
        point=point-2*t*z*M::square(alpha)*(M::square(proj[0])+M::square(proj[1])+proj[0]*proj[1])/(3*M::square(mass));
        point=point-2*t*z*alpha*(h+k)*(proj[0]+proj[1])/(A*mass);
        auto bubble=disjoint_interval::bubble(fourier_e2_root::Z(t));
        if(!arb_contains_zero(acb_imagref(bubble.v)))throw std::runtime_error("real controlled bubble");
        arb_set(B.v,acb_realref(bubble.v));full=M::density(r,chart,p).value;
    }else{
        // Infer the epsilon^2 and epsilon coefficients at original v nodes.
        // Simpson integration is exact for C(v), L(v), and L(v)^2.
        Real cb,lb,l2b;
        for(unsigned j=0;j<3;++j){
            Real v=Real(long(j))/2,weight=Real(j==1?4:1)/6;
            Real plus=O::original_variance(r,t,y,q1,q2,v,1),minus=O::original_variance(r,t,y,q1,q2,v,-1);
            Real C=(plus+minus)/2-mass,L=(plus-minus)/2;
            cb=cb+weight*C;lb=lb+weight*L;l2b=l2b+weight*O::sq(L);
        }
        Real l1=2*(r.host[0][0]*q1+r.host[0][1]*q2),l2=2*(r.host[1][0]*q1+r.host[1][1]*q2);
        point=((O::sq(l1)+O::sq(l2)+l1*l2)/O::sq(A)-2/A)*O::lg(mass)+cb/mass-l2b/(2*O::sq(mass))-(l1+l2)*lb/(A*mass);
        B=O::original_bubble(t);full=O::density(r,chart,p);
    }
    Real common=2/(s*(1-s))*(1-y)*jac*M::power(t/A,3)*(B/A);
    Real result=full*2/(r.ln2*s*(1-s))+common*(1-z)*point;
    if(!arb_is_finite(result.v))throw std::runtime_error("finite full residual required");return result;
}
}
int main(int argc,char**argv){try{using namespace disjoint_control_residual;
    if(argc!=2)throw std::invalid_argument("pinned graph recipe required");
    std::string mode;unsigned bits,order,count;
    if(!(std::cin>>mode>>bits>>order>>count)||(mode!="stable"&&mode!="original")||bits<512||bits>2048||
       (order!=1&&order!=4&&order!=8&&order!=16)||count<1||count>512)throw std::invalid_argument("bounded residual batch required");
    precision=bits;Recipe r;disjoint_interval::prepare(r,argv[1]);O::validate(r);
    std::cout<<"{\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"orbit\":"<<order<<",\"count\":"<<count<<",\"points\":[";
    for(unsigned i=0;i<count;++i){
        std::array<Real,4>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open cube required");}
        if(i)std::cout<<',';Real total;std::cout<<"{\"charts\":[";
        for(unsigned chart=0;chart<2;++chart){
            Real sum;
            for(unsigned j=0;j<order;++j){auto pp=p;pp[2]=p[2]+Real(long(j))/Real(long(order));
                if(arb_ge(pp[2].v,Real(1).v))pp[2]=pp[2]-1;
                else if(!arb_lt(pp[2].v,Real(1).v))throw std::runtime_error("uncertain angular orbit boundary");
                sum=sum+value(r,mode,chart,pp)/Real(long(order));
            }
            if(chart)std::cout<<',';interval(sum);total=total+sum;
        }
        std::cout<<"],\"total\":";interval(total);std::cout<<'}';
    }
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing residual request");
    std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
