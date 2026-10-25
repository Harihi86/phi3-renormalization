// Isolated interval evaluation of original L3E3G6 and its exact scale reduction.
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
struct Chart{Real q,u,n,d,jac7,measure8,ca,cb,x0,x1,original_jac,mass;std::array<Real,9>alpha;};
Chart chart(const std::vector<Real>&p,const std::array<int,3>&perm){
    const Real&s=p[0],&v=p[1],&f0=p[2],&f1=p[3],&f2=p[4],&r0=p[5],&t=p[6];
    Real den=1+s+s*v;std::array<Real,3>w,f{f0,f1,f2};
    std::array<Real,3>ordered{1/den,s/den,s*v/den};for(int i=0;i<3;++i)w[perm[i]]=ordered[i];
    Real psi=w[0]*w[1]+w[0]*w[2]+w[1]*w[2];
    Real q=0;for(int i=0;i<3;++i)q+=w[i]*f[i]*(1-f[i]);q+=w[0]*w[1]*w[2]*norm(f0-f2,f1-f2)/psi;
    Real mu0=w[1]*w[2]*(f0-f2)/psi-f0,mu1=w[1]*w[2]*(f1-f2)/psi;
    Real rd=square(r0)+square(1-r0),r=square(r0)/rd,rc=square(1-r0)/rd,dr=2*r0*(1-r0)/square(rd);
    Chart c;c.q=q;c.u=r*rc*(w[1]+w[2])/psi;c.n=r*rc*norm(mu0+t,mu1);c.d=rc*t*(1-t);
    c.jac7=v/integer_power(1+v+s*v,3)*rc*dr;
    if(p.size()==7)return c;
    Real z0=p[7],z=square(z0),y=(1-z0)*(1+z0),cc=1+c.u,sigma=y/(cc*z);
    Real d1=c.u+z,d2=1+c.u*z,d3=1+z*(c.u+q*cc);
    c.x0=q*cc*z/d2;c.x1=y*(c.n*cc*z/d1+c.d)/d3;
    c.ca=c.u*(square(d1)+d1*cc*z+square(cc*z))/integer_power(d1,3)*q*cc/d2;
    c.cb=integer_power(cc,3)*square(z)/integer_power(d1,3)*(c.n*cc*z/d1+c.d)/d3;
    c.measure8=2*z0*c.jac7;
    c.original_jac=w[0]*w[1]*w[2]*s/integer_power(den,3)*rc*dr*2*z0*square(sigma)/(cc*square(z));
    c.mass=1+sigma;
    for(int i=0;i<3;++i){c.alpha[i]=w[i]*(1-f[i]);c.alpha[i+3]=w[i]*f[i];}
    c.alpha[6]=sigma*r;c.alpha[7]=sigma*rc*t;c.alpha[8]=sigma*rc*(1-t);
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
int main(){try{
    std::string mode;unsigned bits,count;
    if(!(std::cin>>mode>>bits>>count)||bits<128||bits>8192||!count||count>4096)throw std::invalid_argument("bounded mode, precision and batch required");
    if(mode!="K"&&mode!="scale"&&mode!="scale7"&&mode!="stable8"&&mode!="original8")throw std::invalid_argument("unknown mode");
    precision=bits;std::vector<State>states;
    if(mode=="original8")for(unsigned j=0;j<2;++j){int sign;if(!(std::cin>>sign)||(sign!=1&&sign!=-1))throw std::invalid_argument("signed original forest state required");states.push_back({sign,read_poly(),read_poly()});}
    std::array<std::array<int,3>,6>perms{{{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}}};
    std::cout<<"{\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"points\":[";
    for(unsigned index=0;index<count;++index){
        if(index)std::cout<<',';
        std::vector<Real>p;unsigned dim=mode=="K"?1:(mode=="scale"?4:(mode=="scale7"?7:8));
        for(unsigned j=0;j<dim;++j){p.push_back(number());if(dim>=7&&(!arb_is_positive(p.back().v)||!arb_lt(p.back().v,Real(1).v)))throw std::invalid_argument("strict open cube required");}
        if(mode=="K"){interval(K(p[0]));continue;}
        if(mode=="scale"){interval(scale(p[0],p[1],p[2],p[3]));continue;}
        Real total=0;std::cout<<"{\"sectors\":[";
        for(unsigned sector=0;sector<6;++sector){
            Chart c=chart(p,perms[sector]);Real value;
            if(mode=="scale7")value=c.jac7*scale(c.q,c.u,c.n,c.d);
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
