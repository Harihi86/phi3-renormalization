// Generic executor for the serialized recursive Gaussian BPH program.
// All routes, signs, Taylor orders and monomials come from the input program.
// Boost is a separately supplied header-only build dependency; no topology table.
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {
using F50 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50>, boost::multiprecision::et_off>;
using F100 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off>;
thread_local std::string last_error;

struct Layout {
    std::vector<int> orders, strides;
    std::vector<std::vector<std::pair<int,int>>> pairs;
    int size=1, degree=0;
    explicit Layout(std::vector<int> o): orders(std::move(o)) {
        strides.resize(orders.size());
        for(int j=int(orders.size())-1;j>=0;--j) {
            strides[j]=size; size*=orders[j]+1; degree+=orders[j];
        }
        if(size>256) throw std::runtime_error("Taylor layout exceeds native capacity");
        std::vector<std::vector<int>> powers(size,std::vector<int>(orders.size()));
        for(int a=0;a<size;++a) for(size_t j=0;j<orders.size();++j)
            powers[a][j]=(a/strides[j])%(orders[j]+1);
        pairs.resize(size);
        for(int k=0;k<size;++k) for(int a=0;a<size;++a) {
            bool valid=true;int b=0;
            for(size_t j=0;j<orders.size();++j) {
                if(powers[a][j]>powers[k][j]) {valid=false;break;}
                b+=(powers[k][j]-powers[a][j])*strides[j];
            }
            if(valid) pairs[k].push_back({a,b});
        }
    }
};
struct Entry {int edge, index, coefficient;};
struct Term {
    int sign;
    Layout layout;
    std::vector<std::vector<Entry>> matrix_entries;
    Term(int s,std::vector<int> o,int width):sign(s),layout(std::move(o)),matrix_entries(width*width) {}
};
struct Program {int loops,edges,legs,width;std::vector<Term> terms;};

template<class T> struct Jet {
    const Layout* l;
    std::vector<T> c;
    explicit Jet(const Layout* layout,T value=0):l(layout),c(layout->size) {c[0]=value;}
    Jet operator-() const {Jet r(l);for(int i=0;i<l->size;++i)r.c[i]=-c[i];return r;}
    Jet operator+(const Jet& b) const {Jet r(l);for(int i=0;i<l->size;++i)r.c[i]=c[i]+b.c[i];return r;}
    Jet operator-(const Jet& b) const {Jet r(l);for(int i=0;i<l->size;++i)r.c[i]=c[i]-b.c[i];return r;}
    Jet operator*(const T& b) const {Jet r(l);for(int i=0;i<l->size;++i)if(c[i]!=0)r.c[i]=c[i]*b;return r;}
    Jet operator/(const T& b) const {Jet r(l);for(int i=0;i<l->size;++i)if(c[i]!=0)r.c[i]=c[i]/b;return r;}
    Jet operator*(const Jet& b) const {
        Jet r(l);
        for(int k=0;k<l->size;++k) for(auto ab:l->pairs[k])
            if(c[ab.first]!=0 && b.c[ab.second]!=0)r.c[k]+=c[ab.first]*b.c[ab.second];
        return r;
    }
    Jet inverse() const {
        if(c[0]==0)throw std::runtime_error("Singular constant Gaussian pivot");
        // Coefficient convolution a*b=1 is triangular in the mixed-radix
        // monomial index. This avoids repeated full geometric-series products.
        const T inv=T(1)/c[0];Jet result(l,inv);
        for(int k=1;k<l->size;++k) {
            for(auto ab:l->pairs[k])if(ab.first && c[ab.first]!=0 && result.c[ab.second]!=0)
                result.c[k]-=c[ab.first]*result.c[ab.second];
            result.c[k]*=inv;
        }
        return result;
    }
};

template<class T> T stable_log1p(const T& z) {
    using boost::multiprecision::abs;using boost::multiprecision::log;
    if(abs(z)>=T("0.125"))return log(T(1)+z);
    T term=z,result=z;
    // No fixed truncation: the active precision determines stabilization.
    for(int n=2;n<100000;++n) {
        term=term*(-z)*T(n-1)/T(n);T next=result+term;
        if(next==result)return result;result=next;
    }
    throw std::runtime_error("log1p series did not stabilize");
}
template<class T> Jet<T> radial(const Jet<T>& h,const T& mass,int legs) {
    using boost::multiprecision::abs;
    const Layout* l=h.l;T h0=h.c[0],z=h0/mass;
    if(T(1)+z<=0)throw std::runtime_error("Invalid logarithm in radial kernel");
    T logarithm=stable_log1p(z);
    Jet<T> delta=h;delta.c[0]=0;
    if(legs==3) {
        Jet<T> normalized=delta/(mass+h0),term(l,1),result(l,-logarithm);
        for(int n=1;n<=l->degree;++n){term=term*normalized;result=result+term*(T(n%2?-1:1)/T(n));}
        return result;
    }
    T constant;
    if(abs(z)<T("0.125")) {
        T term=mass*z*z/2;constant=term;
        for(int n=3;n<100000;++n) {
            term=term*(-z)*T(n-2)/T(n);T next=constant+term;
            if(next==constant)break;
            constant=next;
            if(n==99999)throw std::runtime_error("Phi2 series did not stabilize");
        }
    } else constant=(mass+h0)*logarithm-h0;
    Jet<T> result(l,constant);result=result+delta*logarithm;
    Jet<T> term=delta;
    T inverse_scale=T(1)/(mass+h0),factor=inverse_scale;
    for(int n=2;n<=l->degree;++n) {
        term=term*delta;
        result=result+term*(T(n%2?-1:1)*factor/T(n*(n-1)));
        factor*=inverse_scale;
    }
    return result;
}

template<class T> struct Evaluation {T value=0, absolute_components=0;};
template<class T> Evaluation<T> evaluate(const Program& p,const double* alpha,const double* gram,double mass_double) {
    using boost::multiprecision::abs;
    Evaluation<T> output;if(p.legs==1)return output;
    std::vector<T> aa(p.edges),gg((p.legs-1)*(p.legs-1));
    for(int e=0;e<p.edges;++e)aa[e]=T(alpha[e]);
    for(size_t j=0;j<gg.size();++j)gg[j]=T(gram[j]);
    const T mass(mass_double);
    for(const Term& term:p.terms) {
        const Layout* l=&term.layout;int w=p.width;
        std::vector<Jet<T>> matrix;matrix.reserve(w*w);
        for(int j=0;j<w*w;++j) {
            matrix.emplace_back(l);
            for(const Entry& entry:term.matrix_entries[j])
                matrix.back().c[entry.index]+=aa[entry.edge]*entry.coefficient;
        }
        Jet<T> determinant(l,1);
        for(int k=0;k<p.loops;++k) {
            const Jet<T> pivot=matrix[k*w+k];determinant=determinant*pivot;
            const Jet<T> inverse=pivot.inverse();
            for(int i=k+1;i<w;++i) {
                const Jet<T> factor=matrix[i*w+k]*inverse;
                for(int j=i;j<w;++j) {
                    matrix[i*w+j]=matrix[i*w+j]-factor*matrix[k*w+j];
                    if(i!=j)matrix[j*w+i]=matrix[i*w+j];
                }
            }
        }
        Jet<T> h(l);
        for(int i=0;i<p.legs-1;++i) for(int j=0;j<p.legs-1;++j)
            h=h+matrix[(p.loops+i)*w+p.loops+j]*gg[i*(p.legs-1)+j];
        Jet<T> inverse=determinant.inverse();
        Jet<T> value=(inverse*inverse*inverse)*radial(h,mass,p.legs);
        for(const T& coefficient:value.c) {
            output.value+=coefficient*term.sign;
            output.absolute_components+=abs(coefficient);
        }
    }
    return output;
}

void write_error(const std::exception& e){last_error=e.what();}
} // namespace

extern "C" {
int phi3lab_bph_abi_version(){return 1;}
const char* phi3lab_bph_last_error(){return last_error.c_str();}
void* phi3lab_bph_create(const int* input,size_t count) {
    try {
        size_t cursor=0;
        auto read=[&](){if(cursor>=count)throw std::runtime_error("Truncated native program");return input[cursor++];};
        if(read()!=1)throw std::runtime_error("Native program ABI mismatch");
        auto p=std::make_unique<Program>();p->loops=read();p->edges=read();p->legs=read();int nt=read();
        if(p->loops<1||p->loops>4||p->edges<1||p->edges>20||p->legs<1||p->legs>3||nt<1||nt>4096)
            throw std::runtime_error("Invalid native program dimensions");
        p->width=p->loops+p->legs-1;
        for(int ti=0;ti<nt;++ti) {
            int sign=read(),nv=read();if((sign!=1&&sign!=-1)||nv<0||nv>6)throw std::runtime_error("Invalid native Taylor term");
            std::vector<int> orders(nv);for(int& o:orders){o=read();if(o<0||o>8)throw std::runtime_error("Invalid Taylor order");}
            p->terms.emplace_back(sign,orders,p->width);Term& term=p->terms.back();
            std::vector<int> coefficients(p->edges*p->width),powers(p->edges*p->width*nv);
            for(int k=0;k<p->edges*p->width;++k) {
                coefficients[k]=read();
                for(int j=0;j<nv;++j){powers[k*nv+j]=read();if(powers[k*nv+j]<0)throw std::runtime_error("Negative monomial power");}
            }
            for(int e=0;e<p->edges;++e)for(int i=0;i<p->width;++i)for(int j=i;j<p->width;++j) {
                int ci=coefficients[e*p->width+i],cj=coefficients[e*p->width+j];if(!ci||!cj)continue;
                int index=0;bool retain=true;
                for(int v=0;v<nv;++v){int power=powers[(e*p->width+i)*nv+v]+powers[(e*p->width+j)*nv+v];if(power>orders[v]){retain=false;break;}index+=power*term.layout.strides[v];}
                if(retain) {
                    Entry entry{e,index,ci*cj};term.matrix_entries[i*p->width+j].push_back(entry);
                    if(i!=j)term.matrix_entries[j*p->width+i].push_back(entry);
                }
            }
        }
        if(cursor!=count)throw std::runtime_error("Trailing native program data");
        return p.release();
    }catch(const std::exception& e){write_error(e);return nullptr;}
}
void phi3lab_bph_destroy(void* handle){delete static_cast<Program*>(handle);}

// status 0: 50/100-digit agreement with a conservative cancellation guard.
// status 1: precision insufficient, 2: arithmetic failure; Python must retry.
// The agreement criterion is numerical evidence, not a rigorous interval bound.
int phi3lab_bph_batch(void* handle,size_t count,const double* alpha,const double* gram,double mass,
                     double* output,int* status,double* relative_check) {
    if(!handle)return -1;
    const Program& p=*static_cast<Program*>(handle);
    for(size_t i=0;i<count;++i) {
        output[i]=std::numeric_limits<double>::quiet_NaN();status[i]=2;relative_check[i]=std::numeric_limits<double>::infinity();
        try {
            auto low=evaluate<F50>(p,alpha+i*p.edges,gram,mass);
            auto high=evaluate<F100>(p,alpha+i*p.edges,gram,mass);
            using boost::multiprecision::abs;
            // Relative comparison protects genuinely tiny nonzero densities.
            // A nontrivial exact cancellation to zero is conservatively retried.
            F100 scale=abs(high.value);
            F100 difference=(scale==0 ? F100(low.value==0?0:1) : F100(abs(high.value-F100(low.value))/scale));
            F100 guard=(scale==0 ? F100(high.absolute_components==0?0:1) : F100(high.absolute_components*F100("1e-90")/scale));
            relative_check[i]=std::max(difference,guard).convert_to<double>();
            output[i]=high.value.convert_to<double>();
            status[i]=(difference<=F100("1e-25")&&guard<=F100("1e-25")&&std::isfinite(output[i]))?0:1;
            if(output[i]==0 && high.value!=0)status[i]=3; // binary64 underflow, no relative guarantee
        }catch(const std::exception& e){write_error(e);}
    }
    return 0;
}
int phi3lab_bph_batch_parallel(void* handle,size_t count,const double* alpha,const double* gram,double mass,
                              int threads,double* output,int* status,double* relative_check) {
    if(!handle)return -1;
    if(count==0)return 0;
    threads=std::max(1,std::min(threads,64));
    threads=static_cast<int>(std::min(count,static_cast<size_t>(threads)));
    if(threads==1)return phi3lab_bph_batch(handle,count,alpha,gram,mass,output,status,relative_check);
    const Program& p=*static_cast<Program*>(handle);
    std::vector<std::thread> workers;workers.reserve(threads);
    try {
        for(int thread=0;thread<threads;++thread) {
            size_t begin=count*thread/threads,end=count*(thread+1)/threads;
            workers.emplace_back([=,&p](){phi3lab_bph_batch(handle,end-begin,alpha+begin*p.edges,gram,mass,
                                                         output+begin,status+begin,relative_check+begin);});
        }
    }catch(const std::exception& e) {
        for(auto& worker:workers)if(worker.joinable())worker.join();write_error(e);return -1;
    }
    for(auto& worker:workers)worker.join();return 0;
}
int phi3lab_bph_text(void* handle,const double* alpha,const double* gram,double mass,int digits,char* output,size_t capacity) {
    try {
        if(!handle)throw std::runtime_error("Null native program");
        const Program& p=*static_cast<Program*>(handle);std::ostringstream stream;
        if(digits==50)stream<<std::setprecision(50)<<evaluate<F50>(p,alpha,gram,mass).value;
        else if(digits==100)stream<<std::setprecision(100)<<evaluate<F100>(p,alpha,gram,mass).value;
        else throw std::runtime_error("Use 50 or 100 native digits");
        std::string value=stream.str();if(value.size()+1>capacity)throw std::runtime_error("Native text output buffer too small");
        std::memcpy(output,value.c_str(),value.size()+1);return 0;
    }catch(const std::exception& e){write_error(e);return -1;}
}
}
