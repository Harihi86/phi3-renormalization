// Three exact scale sectors and analytic integration of one simplex coordinate.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
namespace overlap {
using arb_cw::Real;using arb_cw::precision;
using fourier_e2_root::readq;using fourier_e2_root::interval;
namespace Pow=disjoint_v_integrated;
Real sq(const Real&x){Real r;arb_sqr(r.v,x.v,precision);return r;}
Real rt(const Real&x){Real r;arb_sqrt(r.v,x.v,precision);return r;}
Real log(const Real&x){Real r;arb_log(r.v,x.v,precision);return r;}
Real lp(const Real&x){Real r;arb_log1p(r.v,x.v,precision);return r;}
Real ah(const Real&x){Real r;arb_atanh(r.v,x.v,precision);return r;}
void pos(const Real&x,const char*s){if(!arb_is_positive(x.v))throw std::runtime_error(s);}
Real h(const Real&x){return (1+x)*lp(x)-x;}
Real mean_h(const Real&t,const Real&M,const Real&A,const Real&B,const Real&C,bool constant=false){
    if(constant)return h(t*C/M);
    pos(A,"positive conditional curvature");Real centre=B/(2*A),D=t*A/M,S=1+t*(C+sq(B)/(4*A))/M;
    pos(D,"positive conditional scale");pos(S,"positive conditional maximum");
    Real rootD=rt(D),rootS=rt(S);
    auto primitive=[&](const Real&v){
        Real f=S-D*sq(v),arg=v*rootD/rootS;pos(f,"conditional logarithm domain");pos(1-sq(arg),"conditional atanh domain");
        return (S*v-D*Pow::power(v,3)/3)*log(f)+5*D*Pow::power(v,3)/9+(1-7*S/3)*v+4*S*rootS/(3*rootD)*ah(arg);
    };
    return primitive(1-centre)-primitive(-centre);
}
struct UF{Real u,f;};
UF original_uf(const Real&a,const Real&b,const Real&x,const Real&y,unsigned state){
    // Frozen original Gaussian routes, in the common affine loop basis (k,l,q).
    static const int routes[3][5][3]={
        {{1,1,1},{-1,-1,0},{-1,0,-1},{1,0,0},{0,1,0}},
        {{0,1,0},{0,-1,0},{-1,0,-1},{1,0,0},{0,1,0}},
        {{-1,0,1},{1,0,0},{0,1,0},{0,-1,0},{0,1,0}}};
    std::array<Real,5> alpha={a*x,a*(1-x),b*y,b*(1-y),Real(1)};
    Real m[3][3];
    for(unsigned e=0;e<5;++e)for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
        m[i][j]=m[i][j]+alpha[e]*Real(long(routes[state][e][i]*routes[state][e][j]));
    Real u=m[0][0]*m[1][1]-m[0][1]*m[1][0];
    Real f=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
           -m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
           +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
    pos(u,"positive original Gaussian determinant");return {u,f};
}
Real core(const std::string&mode,const Real&t,const Real&a,const Real&b,const Real&y){
    Real M=1+a+b;
    if(mode=="closed"){
        Real U=a*b+a+b,A=sq(a)*(b+1)/U,B=a-2*a*b*y/U,C=b*y*(1-y)+a*b*sq(y)/U;
        return a*b*M*(mean_h(t,M,A,B,C)/Pow::power(U,3)
            -mean_h(t,M,Real(0),Real(0),b*y*(1-y),true)/Pow::power(b*(a+1),3)
            -mean_h(t,M,a,a,Real(0))/Pow::power(a*(b+1),3));
    }
    Real sum;
    for(unsigned state=0;state<3;++state){
        auto p0=original_uf(a,b,Real(0),y,state),pm=original_uf(a,b,Real(1)/2,y,state),p1=original_uf(a,b,Real(1),y,state);
        if(!arb_overlaps(p0.u.v,pm.u.v)||!arb_overlaps(p0.u.v,p1.u.v))throw std::runtime_error("original U not constant in x");
        Real q0=p0.f/p0.u,qm=pm.f/pm.u,q1=p1.f/p1.u,A=2*(2*qm-q0-q1),B=q1-q0+A;
        if(state==1&&(!arb_contains_zero(A.v)||!arb_contains_zero(B.v)))throw std::runtime_error("original constant forest changed");
        Real one=mean_h(t,M,A,B,q0,state==1)/Pow::power(p0.u,3);
        sum=sum+(state==0?one:-one);
    }
    return a*b*M*sum;
}
Real outer(const std::string&gid,const Real&t,const Real&nu){
    Real A=1+t;
    if(gid=="L3E2G7"){
        Real delta=2/(t+rt(sq(t)+4)),e=A+delta;
        return delta*(3*sq(A)+6*A*delta+(3+t)*sq(delta))/(3*Pow::power(A,3)*Pow::power(e,3));
    }
    Real delta=2*(t*nu*(1-nu)+1)/(rt(sq(t)+4*t*nu*(1-nu)+4)+t);
    return -delta*(2*A+delta)/(sq(A)*sq(A+delta));
}
Real value(const std::string&gid,const std::string&mode,unsigned tp,unsigned control,const std::array<Real,5>&p){
    const Real&s=p[0];const Real&r=p[1];const Real&v=p[2];Real y=p[3]/2,nu=p[4]/2;
    Real t=Pow::power(s/(1-s),tp),A=1+t,dt=Real(long(tp))*t/(s*(1-s));
    Real density=2*r*core(mode,t,r*v,r,y)+2/sq(v)*core(mode,t,r,1/v,y)
                 +2/(Pow::power(r,3)*sq(v))*core(mode,t,1/r,1/(r*v),y);
    if(control){
        // Exact conditional integral of the r->0 third-sector limit.
        Real M=1+v,constant=y*(1-y);
        Real H=mean_h(t,M,v,v,constant)-h(t*constant/M)-mean_h(t,M,v,v,Real(0));
        Real limit=2*(1+v)/v*H;
        density=density+limit*(Real(1)/Real(long(control+1))-Pow::power(1-r,control));
    }
    Real out=sq(t)/(2*sq(A))*dt*outer(gid,t,nu)*density;
    if(!arb_is_finite(out.v))throw std::runtime_error("finite whole overlap density required");return out;
}
std::string recipe(const char*path){
    std::ifstream in(path);std::string id,sha,extra;
    if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("exact graph recipe required");
    bool valid=(id=="L3E2G7"&&sha=="97ed1d0710c961669079e1f5665910f1442147f2ac0034b060201f72d55d1b97")||
      (id=="L3E3G38"&&sha=="a3b16b71c29c32a18e858f3ae257b956ff1b7b5271b136eafe731387e0f9ed50");
    if(!valid)throw std::invalid_argument("frozen overlapping graph required");return id;
}
}
int main(int argc,char**argv){try{using namespace overlap;
    if(argc!=2)throw std::invalid_argument("one graph recipe required");auto gid=recipe(argv[1]);
    std::string mode;unsigned bits,tp,control,count;
    if(!(std::cin>>mode>>bits>>tp>>control>>count)||(mode!="closed"&&mode!="original")||bits<512||bits>2048||
       (tp!=1&&tp!=2)||control>2||count<1||count>512)throw std::invalid_argument("bounded overlap batch required");
    precision=bits;
    std::cout<<"{\"id\":\""<<gid<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"t_power\":"<<tp<<",\"control\":"<<control<<",\"sectors\":3,\"analytic_x\":true,\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){
        std::array<Real,5>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open five-cube required");}
        auto out=value(gid,mode,tp,control,p);if(i)std::cout<<',';interval(out);
    }
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing overlap request");
    std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
