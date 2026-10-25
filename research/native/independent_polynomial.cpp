// Independent positive graph-polynomial / squarefree BPH evaluator.
// No production Gaussian, routing, rectangular jet, or counterterm code.
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <vector>
#include <cmath>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <thread>
#include <atomic>
#include <cstdlib>
using R50 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>>;
using R100 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>>;
struct Term { int mask, slot, channel; };
struct State { int k, sign; std::vector<Term> u, f; };
struct Program { int edges, legs; std::vector<State> states; std::vector<int> masks; };
template<class R> using Series = std::vector<R>;
template<class R> Series<R> mul(const Series<R>& a,const Series<R>& b) {
    Series<R> out(a.size(),R(0));
    for(unsigned m=0;m<a.size();++m) {
        unsigned sub=m;
        while(true) { out[m]+=a[sub]*b[m^sub]; if(!sub)break;sub=(sub-1)&m; }
    }
    return out;
}
template<class R> Series<R> inverse(const Series<R>& a,int k) {
    Series<R> h(a.size(),R(0)), out(a.size(),R(0)), power(a.size(),R(0));
    out[0]=power[0]=1;
    for(unsigned m=1;m<a.size();++m)h[m]=a[m]/a[0];
    for(int degree=1;degree<=k;++degree) {
        power=mul(power,h);
        for(unsigned m=0;m<a.size();++m)out[m]+=(degree%2?-power[m]:power[m]);
    }
    for(auto& v:out)v/=a[0];
    return out;
}
template<class R> R phi2(const R& h,const R& mass) {
    R x=h/mass;
    if(abs(x)>=R("0.05"))return (mass+h)*log1p(x)-h;
    R term=x*x,total=term/2;
    const R eps=std::numeric_limits<R>::epsilon()/8;
    for(int n=3;n<2000;++n) {
        term*=-x;R addition=term/(n*(n-1));total+=addition;
        if(addition==0 || abs(addition)<=eps*abs(total))return mass*total;
    }
    throw std::runtime_error("Phi2 series did not converge");
}
template<class R> Series<R> density(const Series<R>& u,const Series<R>& f,int k,int legs,const R& mass) {
    Series<R> inv=inverse(u,k),h=mul(f,inv),delta=h,power(u.size(),R(0)),phi(u.size(),R(0));
    delta[0]=0;power[0]=1;
    std::vector<R> coeff(k+1);
    if(legs==3) {
        coeff[0]=-log1p(h[0]/mass);
        for(int n=1;n<=k;++n)coeff[n]=(n%2?R(-1):R(1))/(R(n)*pow(mass+h[0],n));
    } else {
        coeff[0]=phi2(h[0],mass);
        if(k)coeff[1]=log1p(h[0]/mass);
        for(int n=2;n<=k;++n)coeff[n]=(n%2?R(-1):R(1))/(R(n*(n-1))*pow(mass+h[0],n-1));
    }
    phi[0]=coeff[0];
    for(int n=1;n<=k;++n) {
        power=mul(power,delta);
        for(unsigned m=0;m<phi.size();++m)phi[m]+=coeff[n]*power[m];
    }
    return mul(mul(mul(inv,inv),inv),phi);
}
template<class R> std::pair<R,R> evaluate(const Program& p,const double* alpha,const double* invariants,double mass) {
    std::vector<R> monomial(1<<p.edges,R(0));
    for(int mask:p.masks) {
        R product=1;
        for(int edge=0;edge<p.edges;++edge)if(mask&(1<<edge))product*=R(alpha[edge]);
        monomial[mask]=product;
    }
    R total=0,absolute=0;
    for(const auto& state:p.states) {
        Series<R> u(1<<state.k,R(0)),f(1<<state.k,R(0));
        for(const auto& t:state.u)u[t.slot]+=monomial[t.mask];
        for(const auto& t:state.f)f[t.slot]+=monomial[t.mask]*R(invariants[t.channel]);
        auto v=density(u,f,state.k,p.legs,R(mass));
        for(const R& x:v) {total+=state.sign*x;absolute+=abs(x);}
    }
    return {total,absolute};
}
extern "C" {
int independent_polynomial_abi(){return 1;}
void* independent_polynomial_create(const int* data,int length) {
    try {
        int j=0;
        auto next=[&](){if(j>=length)throw std::runtime_error("truncated program");return data[j++];};
        auto p=new Program();p->edges=next();p->legs=next();int count=next();
        if(p->edges<1 || p->edges>18 || (p->legs!=2&&p->legs!=3))throw std::runtime_error("invalid graph scope");
        std::vector<int> masks;
        for(int state_index=0;state_index<count;++state_index) {
            State s;s.k=next();s.sign=next();int nu=next(),nf=next();
            if(s.k<0||s.k>8)throw std::runtime_error("unsupported Taylor size");
            for(int t=0;t<nu+nf;++t) {
                Term term{next(),next(),next()};
                if(term.mask<0||term.mask>=(1<<p->edges)||term.slot<0||term.slot>=(1<<s.k))throw std::runtime_error("invalid term");
                if(t<nu)s.u.push_back(term);else s.f.push_back(term);
                masks.push_back(term.mask);
            }
            p->states.push_back(s);
        }
        if(j!=length)throw std::runtime_error("program trailing integers");
        std::sort(masks.begin(),masks.end());masks.erase(std::unique(masks.begin(),masks.end()),masks.end());p->masks=masks;
        return p;
    } catch(...) {return nullptr;}
}
void independent_polynomial_free(void* handle){delete static_cast<Program*>(handle);}
int independent_polynomial_batch(void* handle,int count,const double* alpha,const double* invariants,double mass,double* output,int* status,double* condition) {
    if(!handle)return -1;
    const auto& p=*static_cast<Program*>(handle);
    std::atomic<int> next(0);
    auto worker=[&](){while(true){
      int start=next.fetch_add(16);if(start>=count)break;
      for(int row=start;row<std::min(count,start+16);++row) {
        try {
            auto a=evaluate<R50>(p,alpha+row*p.edges,invariants,mass);
            // Fifty decimal digits are ample unless the recursive sum loses
            // at least twenty digits. Recompute only those rows at100 digits.
            R50 threshold=R50("1e-25")*abs(a.first);
            if(R50("1e-44")*a.second<=threshold) {
                output[row]=a.first.convert_to<double>();status[row]=50;
            } else {
                auto b=evaluate<R100>(p,alpha+row*p.edges,invariants,mass);
                R100 old(a.first.str(0,std::ios_base::scientific));
                R100 scale=std::max(R100(abs(b.first)),R100(abs(old)));
                if((scale==0&&b.second==0)||(scale!=0&&abs(b.first-old)<=R100("1e-25")*scale&&R100("1e-94")*b.second<=R100("1e-25")*abs(b.first))) {
                    output[row]=b.first.convert_to<double>();status[row]=100;
                } else {output[row]=0;status[row]=-100;}
            }
            condition[row]=(a.first==0?std::numeric_limits<double>::infinity():(a.second/abs(a.first)).template convert_to<double>());
            if(!std::isfinite(output[row]))status[row]=-200;
        } catch(...) {output[row]=0;condition[row]=std::numeric_limits<double>::infinity();status[row]=-300;}
      }
    }};
    const char* requested=std::getenv("PHI3LAB_INDEPENDENT_THREADS");
    int threads=requested?std::atoi(requested):4;
    threads=std::max(1,std::min(threads,std::max(1,(count+31)/32)));
    std::vector<std::thread> pool;
    for(int j=1;j<threads;++j)pool.emplace_back(worker);
    worker();for(auto& thread:pool)thread.join();
    return 0;
}
}
