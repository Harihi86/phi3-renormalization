// Isolated interval evaluation of original L3E3G35 and its exact scale reduction.
#include "arb-cw-native-v2/native_arb.hpp"
#include "fmpq.h"
#include <array>
#include <iostream>
using arb_cw::Real;
using arb_cw::precision;
using arb_cw::integer_power;

Real square(const Real& a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real norm(const Real& x,const Real& y){return square(x-y/2)+3*square(y)/4;}
Real number(){
    std::string n,d;
    if(!(std::cin>>n>>d)||n.size()>10000||d.size()>10000)throw std::invalid_argument("bounded rational required");
    fmpq_t q;fmpq_init(q);
    if(fmpq_set_str(q,(n+"/"+d).c_str(),10)||fmpz_sgn(fmpq_denref(q))<=0){fmpq_clear(q);throw std::invalid_argument("canonical positive denominator required");}
    fmpq_canonicalise(q);Real r;arb_set_fmpq(r.v,q,precision);fmpq_clear(q);return r;
}
void integer(const fmpz_t z){char*s=fmpz_get_str(nullptr,10,z);std::cout<<'"'<<s<<'"';flint_free(s);}
void interval(const Real&r){
    if(!r.finite())throw std::domain_error("nonfinite point enclosure");
    fmpz_t lo,hi,ex;fmpz_init(lo);fmpz_init(hi);fmpz_init(ex);
    arb_get_interval_fmpz_2exp(lo,hi,ex,r.v);
    std::cout<<'[';integer(lo);std::cout<<',';integer(hi);std::cout<<',';integer(ex);std::cout<<']';
    fmpz_clear(lo);fmpz_clear(hi);fmpz_clear(ex);
}
void positive(const Real&r){if(!arb_is_positive(r.v))throw std::domain_error("strictly positive input required");}
void nonnegative(const Real&r){if(!arb_is_nonnegative(r.v))throw std::domain_error("nonnegative input required");}
Real log_ratio(const Real&x){return x.zero()?Real(1):arb_cw::log1p(x)/x;}
Real K(const Real&v){
    positive(v);Real pi;arb_const_pi(pi.v,precision);Real constant=square(pi)/6-Real(5)/4;
    if(arb_is_one(v.v))return constant;
    Real delta=v-1,upper=arb_cw::abs_upper(delta);
    if(arb_le(upper.v,(Real(1)/2).v)){
        unsigned terms=unsigned(precision)+8;
        if(!upper.zero()){
            slong exponent=arf_abs_bound_lt_2exp_si(arb_midref(upper.v));
            slong decay=std::max(slong(1),-exponent);
            terms=std::max(24u,unsigned((precision+decay-1)/decay)+8);
        }
        Real result=constant,power=delta;
        for(unsigned j=1;j<=terms;++j){
            result+=Real(j%2?2:-2)*power/(Real(long(j))*Real(long(j))*Real(long(j+1))*Real(long(j+2)));
            power*=delta;
        }
        Real tail=2*integer_power(upper,terms+1)/(Real(long(terms+1))*Real(long(terms+1))*Real(long(terms+2))*Real(long(terms+3))*(1-upper));
        arb_cw::add_error(result,arb_cw::abs_upper(tail));return result;
    }
    Real dilog,arg=1-v;arb_polylog_si(dilog.v,2,arg.v,precision);
    return square(pi)/6-dilog-v/(2*(1-v))+v*(2-3*v)*arb_cw::log(v)/(2*square(1-v));
}
Real scale(const Real&q,const Real&u,const Real&n,const Real&d){
    nonnegative(q);positive(u);nonnegative(n);nonnegative(d);
    Real a=1+q,x=(1+d)/(a*u),h=n/(a*u);
    Real disc=square(x-1)+2*h*(1+x)+square(h);
    Real plus=(1+x+h+arb_cw::sqrt(disc))/2,minus=x/plus,la=arb_cw::log1p(q);
    return square(la)/2+(arb_cw::log(u)+Real(3)/2)*la-K(plus)-K(minus)+K(1/u)+K(Real(1));
}
struct Chart{Real q,u,n,d,jac7,projective7,measure8,ca,cb,x0,x1,original_jac,mass;std::array<Real,9>alpha;};
Chart chart(const std::vector<Real>&p,const std::array<int,3>&perm){
    const Real&s=p[0],&v=p[1],&f0=p[2],&f1=p[3],&f2=p[4],&r0=p[5],&t=p[6];
    Real den=1+s+s*v;std::array<Real,3>w,f{f0,f1,f2};
    std::array<Real,3>ordered{1/den,s/den,s*v/den};for(int i=0;i<3;++i)w[perm[i]]=ordered[i];
    Real psi=w[0]*w[1]+w[0]*w[2]+w[1]*w[2];
    std::array<Real,3>h;for(int i=0;i<3;++i)h[i]=w[i]*f[i]*(1-f[i]);
    Real kappa=w[0]*w[1]*w[2]/psi,x=f0-f1,y=f0-f2;
    Real a=h[0]+h[1]+kappa*square(x),cross=h[0]+kappa*x*y;
    Real determinant=h[0]*h[1]+h[0]*h[2]+h[1]*h[2]+kappa*(h[0]*square(y-x)+h[1]*square(y)+h[2]*square(x));
    Real rd=square(r0)+square(1-r0),r=square(r0)/rd,rc=square(1-r0)/rd,dr=2*r0*(1-r0)/square(rd);
    std::array<Real,3>host{r,rc*t,rc*(1-t)};
    Chart c;c.q=host[0]*host[1]+host[0]*host[2]+host[1]*host[2];
    c.u=a;c.n=a*norm(cross/a-host[0]-host[2],cross/a-host[0]);c.d=determinant/a;
    c.jac7=v/integer_power(1+v+s*v,3)*rc*dr;
    c.projective7=w[0]*w[1]*w[2]*s/integer_power(den,3)*rc*dr;
    for(int i=0;i<3;++i){c.alpha[i]=w[i]*f[i];c.alpha[i+3]=w[i]*(1-f[i]);c.alpha[i+6]=host[i];}
    if(p.size()==7)return c;
    Real z0=p[7],z=square(z0),y0=(1-z0)*(1+z0),cc=1+c.u,sigma=y0/(cc*z);
    Real d1=c.u+z,d2=1+c.u*z,d3=1+z*(c.u+c.q*cc);
    c.x0=c.q*cc*z/d2;c.x1=y0*(c.n*cc*z/d1+c.d)/d3;
    c.ca=c.u*(square(d1)+d1*cc*z+square(cc*z))/integer_power(d1,3)*c.q*cc/d2;
    c.cb=integer_power(cc,3)*square(z)/integer_power(d1,3)*(c.n*cc*z/d1+c.d)/d3;
    c.measure8=2*z0*c.jac7;
    c.original_jac=w[0]*w[1]*w[2]*s/integer_power(den,3)*rc*dr*2*z0*integer_power(sigma,5)/(cc*square(z));
    c.mass=1+sigma;
    for(int i=0;i<3;++i){c.alpha[i]=sigma*w[i]*f[i];c.alpha[i+3]=sigma*w[i]*(1-f[i]);c.alpha[i+6]=host[i];}
    return c;
}

struct Monomial{unsigned mask;Real coefficient;};
struct State{int sign;std::vector<Monomial>u,f;};
std::vector<Monomial>read_poly(){
    unsigned count;if(!(std::cin>>count)||!count||count>512)throw std::invalid_argument("bounded nonempty polynomial required");
    std::vector<Monomial>result;
    for(unsigned i=0;i<count;++i){unsigned mask;if(!(std::cin>>mask)||mask>=512)throw std::invalid_argument("original edge mask required");result.push_back({mask,number()});}
    return result;
}
Real polynomial(const std::vector<Monomial>&poly,const std::array<Real,9>&alpha){
    Real sum=0;for(const auto&m:poly){Real term=m.coefficient;for(unsigned i=0;i<9;++i)if(m.mask&(1u<<i))term*=alpha[i];sum+=term;}return sum;
}
// Extract the child-scale power directly from each original edge mask.
// This path does not use the Schur chart's Q,u,N,d or simplified jac7.
std::array<Real,5> scale_coefficients(const std::vector<Monomial>&poly,const std::array<Real,9>&alpha,unsigned lo,unsigned hi){
    std::array<Real,5> result;for(Real&r:result)r=0;
    for(const auto&m:poly){
        unsigned degree=__builtin_popcount(m.mask&63u);
        if(degree<lo||degree>hi)throw std::invalid_argument("unexpected original child-scale degree");
        Real term=m.coefficient;for(unsigned i=0;i<9;++i)if(m.mask&(1u<<i))term*=alpha[i];
        result[degree]+=term;
    }
    return result;
}
Real original_scale(const Chart&c,const std::vector<State>&states){
    if(states.size()!=2||states[0].sign!=1||states[1].sign!=-1)throw std::invalid_argument("original bare/CT ordering required");
    auto u=scale_coefficients(states[0].u,c.alpha,2,3),f=scale_coefficients(states[0].f,c.alpha,2,4);
    auto uc=scale_coefficients(states[1].u,c.alpha,2,2),fc=scale_coefficients(states[1].f,c.alpha,2,2);
    if(!arb_overlaps(u[2].v,uc[2].v)||!arb_overlaps(f[2].v,fc[2].v))throw std::domain_error("original CT scale mismatch");
    positive(u[2]);positive(u[3]);Real q=f[2]/u[2],a=u[3]/u[2],d=f[4]/u[3],n=f[3]/u[2]-a*q-d;
    return c.projective7/integer_power(u[2],3)*scale(q,a,n,d);
}
int main(){try{
    std::string mode;unsigned bits,count;
    if(!(std::cin>>mode>>bits>>count)||bits<128||bits>8192||!count||count>4096)throw std::invalid_argument("bounded mode, precision and batch required");
    if(mode!="K"&&mode!="scale"&&mode!="scale7"&&mode!="original7"&&mode!="stable8"&&mode!="original8")throw std::invalid_argument("unknown mode");
    precision=bits;std::vector<State>states;
    if(mode=="original8"||mode=="original7")for(unsigned j=0;j<2;++j){int sign;if(!(std::cin>>sign)||(sign!=1&&sign!=-1))throw std::invalid_argument("signed original forest state required");states.push_back({sign,read_poly(),read_poly()});}
    std::array<std::array<int,3>,6>perms{{{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}};
    std::cout<<"{\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"points\":[";
    for(unsigned index=0;index<count;++index){
        if(index)std::cout<<',';
        std::vector<Real>p;unsigned dim=mode=="K"?1:(mode=="scale"?4:((mode=="scale7"||mode=="original7")?7:8));
        for(unsigned j=0;j<dim;++j){p.push_back(number());if(dim>=7&&(!arb_is_positive(p.back().v)||!arb_lt(p.back().v,Real(1).v)))throw std::invalid_argument("strict open cube required");}
        if(mode=="K"){interval(K(p[0]));continue;}
        if(mode=="scale"){interval(scale(p[0],p[1],p[2],p[3]));continue;}
        Real total=0;std::cout<<"{\"sectors\":[";
        for(unsigned sector=0;sector<6;++sector){
            Chart c=chart(p,perms[sector]);Real value;
            if(mode=="scale7")value=c.jac7*scale(c.q,c.u,c.n,c.d);
            else if(mode=="original7")value=original_scale(c,states);
            else if(mode=="stable8")value=c.measure8*(c.ca*log_ratio(c.x0)-c.cb*log_ratio(c.x1));
            else{
                value=0;for(const auto&state:states){Real u=polynomial(state.u,c.alpha),f=polynomial(state.f,c.alpha);positive(u);nonnegative(f);value+=Real(-state.sign)*arb_cw::log1p(f/u/c.mass)/integer_power(u,3);}
                value*=c.original_jac;
            }
            if(sector)std::cout<<',';interval(value);total+=value;
        }
        std::cout<<"],\"sum\":";interval(total);std::cout<<'}';
    }
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing request input");
    std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
