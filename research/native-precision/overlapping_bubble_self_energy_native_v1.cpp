// Complete weighted-kite forest, six scale charts and exact conditional control.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace weighted_overlap {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
namespace Pow=disjoint_v_integrated;
Real sq(const Real&x){Real r;arb_sqr(r.v,x.v,precision);return r;}
Real rt(const Real&x){Real r;arb_sqrt(r.v,x.v,precision);return r;}
Real ah(const Real&x){Real r;arb_atanh(r.v,x.v,precision);return r;}
void pos(const Real&x,const char*s){if(!arb_is_positive(x.v))throw std::runtime_error(s);}
Real bubble(const Real&t,bool original){
    pos(t,"strictly positive bubble argument");
    if(original){auto value=disjoint_interval::bubble(fourier_e2_root::Z(t));
        if(!arb_contains_zero(acb_imagref(value.v)))throw std::runtime_error("real bubble required");
        Real result;arb_set(result.v,acb_realref(value.v));return result;}
    Real r=rt(t/(t+4));pos(1-r,"bubble atanh interior");
    return (t+4)/3*(ah(r)/r-1)-t/9;
}
struct UF{Real u,q;};
UF original_uf(const Real&a,const Real&b,const Real&x,const Real&y,unsigned state){
    // Effective routes read from the G5 bubble Schur complement: (k,v,p).
    static const int routes[3][5][3]={
        {{0,1,-1},{0,-1,0},{-1,0,1},{1,0,0},{1,-1,0}},
        {{0,1,0},{0,-1,0},{-1,0,1},{1,0,0},{0,-1,0}},
        {{0,1,-1},{0,-1,0},{-1,0,0},{1,0,0},{1,0,0}}};
    std::array<Real,5> alpha={a*x,a*(1-x),b*y,b*(1-y),Real(1)};Real m[3][3];
    for(unsigned e=0;e<5;++e)for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
        m[i][j]=m[i][j]+alpha[e]*Real(long(routes[state][e][i]*routes[state][e][j]));
    Real u=m[0][0]*m[1][1]-m[0][1]*m[1][0];
    Real f=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
    pos(u,"positive original determinant");return {u,f/u};
}
Real density(const std::string&mode,const Real&a,const Real&b,const Real&x,const Real&y,const Real&u,bool control){
    bool original=mode=="original";Real M=1+a+b,d=b*(1-y);pos(M-d,"positive remaining mass");
    std::array<Real,3> U,Q,Q2mean;
    if(!original){
        U={a*b+a+b,b*(a+1),a*(b+1)};
        Real A=sq(a)*(b+1)/U[0],B=a-2*a*b*y/U[0],C=b*y*(1-y)+a*b*sq(y)/U[0];
        Q={-A*sq(x)+B*x+C,b*y*(1-y),a*x*(1-x)};
        Q2mean={sq(A)/5-A*B/2+(sq(B)-2*A*C)/3+B*C+sq(C),sq(Q[1]),sq(a)/30};
    }else{
        for(unsigned state=0;state<3;++state){
            auto point=original_uf(a,b,x,y,state);U[state]=point.u;Q[state]=point.q;
            auto p0=original_uf(a,b,Real(0),y,state),pm=original_uf(a,b,Real(1)/2,y,state),p1=original_uf(a,b,Real(1),y,state);
            Real A=2*(2*pm.q-p0.q-p1.q),B=p1.q-p0.q+A,C=p0.q;
            Q2mean[state]=sq(A)/5-A*B/2+(sq(B)-2*A*C)/3+B*C+sq(C);
        }
    }
    Real sum,point,average;
    for(unsigned state=0;state<3;++state){
        Real mass=M+u*Q[state];pos(mass,"positive auxiliary mass");
        Real j=d*sq(Q[state])*(1-u)*bubble((mass-d)/d,original)/sq(mass),den=Pow::power(U[state],3);
        Real sign=state==0?Real(1):Real(-1);sum=sum+sign*j/den;
        point=point+sign*sq(Q[state])/den;average=average+sign*Q2mean[state]/den;
    }
    if(control){Real f=d*bubble((M-d)/d,original)/sq(M);sum=sum+f*(average/2-(1-u)*point);}
    return a*b*sum;
}
Real value(const std::string&mode,bool control,const std::array<Real,5>&p){
    const auto&r=p[0];const auto&v=p[1];const auto&x=p[2];const auto&y=p[3];const auto&u=p[4];
    std::array<Real,3>aa={r*v,r,1/r},bb={r,1/v,1/(r*v)},jac={r,1/sq(v),1/(Pow::power(r,3)*sq(v))};
    Real total;
    for(unsigned j=0;j<3;++j)for(unsigned orientation=0;orientation<2;++orientation){
        Real a=orientation?bb[j]:aa[j],b=orientation?aa[j]:bb[j];
        total=total+jac[j]*(density(mode,a,b,x,y,u,control)+density(mode,a,b,x,y,1-u,control))/2;
    }
    if(!arb_is_finite(total.v))throw std::runtime_error("finite weighted overlap density required");return total;
}
void recipe(const char*path){
    std::ifstream in(path);std::string id,sha,extra;
    if(!(in>>id>>sha)||id!="L3E2G5"||sha!="9c11f423450cbdd0cb4bb15931d6d9e6cbb8558b11ce31da1f6f270d272335b7"||(in>>extra))
        throw std::invalid_argument("frozen weighted-overlap recipe required");
}
}
int main(int argc,char**argv){try{using namespace weighted_overlap;
    if(argc!=2)throw std::invalid_argument("one graph recipe required");recipe(argv[1]);
    std::string mode;unsigned bits,control,count;
    if(!(std::cin>>mode>>bits>>control>>count)||(mode!="closed"&&mode!="original")||bits<512||bits>2048||control>1||count<1||count>512)
        throw std::invalid_argument("bounded weighted overlap batch required");
    precision=bits;
    std::cout<<"{\"id\":\"L3E2G5\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"control\":"<<control<<",\"sectors\":6,\"auxiliary_antithetic\":2,\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){
        std::array<Real,5>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open five-cube required");}
        auto out=value(mode,control,p);if(i)std::cout<<',';interval(out);
    }
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing weighted overlap request");
    std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
