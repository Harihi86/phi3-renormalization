// Spectral bubble elimination, followed by an analytic common Schwinger radius.
// The first mass derivative is propagated as an Arb jet, never a finite difference.
#include "nested_spectral_scalar.hpp"
namespace nested_spectral {
using namespace disjoint_triangle;
struct Graph {std::string id,sha;std::vector<std::vector<int>> host;std::vector<std::vector<std::vector<int>>> delta;std::vector<std::pair<int,int>> U,F;};
#include "nested_spectral_data.hpp"
struct Jet {Real v,d;Jet():v(),d(){} Jet(long n):v(n),d(){} Jet(const Real&x):v(x),d(){} Jet(const Real&x,const Real&dx):v(x),d(dx){}};
Jet operator+(const Jet&a,const Jet&b){return {a.v+b.v,a.d+b.d};}
Jet operator-(const Jet&a,const Jet&b){return {a.v-b.v,a.d-b.d};}
Jet operator-(const Jet&a){return {-a.v,-a.d};}
Jet operator*(const Jet&a,const Jet&b){return {a.v*b.v,a.d*b.v+a.v*b.d};}
Jet operator/(const Jet&a,const Jet&b){pos(b.v,"positive jet denominator");return {a.v/b.v,(a.d*b.v-a.v*b.d)/sq(b.v)};}
Jet jlog(const Jet&a){pos(a.v,"positive jet logarithm");return {logarithm(a.v),a.d/a.v};}
Jet jroot(const Jet&a){Real r=root(a.v);pos(r,"separated mass-jet roots");return {r,a.d/(2*r)};}
Jet jphi(const Jet&a){auto m=radial_moments(a.v);return {radial_phi(a.v),(m[0]-2*m[1]+m[2])*a.d};}
std::array<Jet,4> jmoments(const Jet&h){
    // Here h<=1/4, so 1-h is uniformly positive.
    Jet z=1-h;pos(z.v,"mass-jet moment recurrence domain");std::array<Jet,4>out;
    out[0]=-jlog(h)/z;for(int n=2;n<=4;++n)out[n-1]=(out[n-2]-Jet(Real(1)/Real(n-1)))/z;return out;
}
struct Coeff {Real q,h,de,ge,mean;};
Jet primitive(const Coeff&c,const Jet&mu,bool control){
    pos(c.h,"positive child variance");Jet h=Jet(c.h)/mu,de=Jet(c.de)/mu,ge=Jet(c.ge)/(mu*mu),mean=Jet(c.mean)/mu;
    Real variance=c.h*c.de-c.ge;
    if(arb_is_negative(variance.v))throw std::runtime_error("inconsistent unscaled Schur variance");
    arb_nonnegative_part(variance.v,variance.v);
    Jet a(1+c.q),b=1+h*a+de,cc=h+ge,delta=1-h*a+de;
    Jet disc=delta*delta+Jet(4*variance)*a/(mu*mu),s=jroot(disc),plus=(b+s)/(2*a),minus=2*cc/(b+s);
    Real la=stable_log1p(c.q);
    Jet answer=Jet(la)*jlog(h)+Jet(sq(la)/2+3*la/2)-jphi(plus/h)-jphi(minus/h)+jphi(Jet(1)/h)+Jet(radial_phi(Real(1)));
    if(control){auto m=jmoments(h);answer=answer+Jet(Real(1)/4-c.q)*h*(m[0]+m[1]+m[2])-(mean-de)*m[3];}
    return answer;
}
std::array<Real,3> polynomial(const std::vector<std::pair<int,int>>&rows,const std::array<Real,3>&alpha,
        const std::array<std::array<Real,3>,2>&child,const std::array<Real,2>&y,int active,bool mean){
    static const int ij[3][2]={{0,1},{0,2},{1,2}};std::array<Real,3>out;
    for(const auto&term:rows){int mask=term.first,degree=0;Real value(term.second);bool use=true;
        if(mean){int d=__builtin_popcount(unsigned(mask&7)),fact=1;for(int j=2;j<=d+2;++j)fact*=j;value=value*2/fact;}
        else for(int e=0;e<3;++e)if(mask&(1<<e))value=value*alpha[e];
        for(int i=0;i<2;++i)for(int j=0;j<3;++j)if(mask&(1<<(3+3*i+j))){
            if(!(active&(1<<i))){use=false;break;}++degree;value=value*y[i]*child[i][ij[j][0]]*child[i][ij[j][1]];}
        if(use){if(degree>2)throw std::runtime_error("invalid reduced graph polynomial degree");out[degree]=out[degree]+value;}
    }return out;
}
std::array<Coeff,3> coefficients(const Graph&g,bool original,const std::array<Real,3>&alpha,
        const std::array<std::array<Real,3>,2>&child,const std::array<Real,2>&y){
    std::array<Coeff,3>out;Real host[3][3];
    for(int e=0;e<3;++e)for(int i=0;i<3;++i)for(int j=0;j<3;++j)host[i][j]=host[i][j]+alpha[e]*Real(g.host[e][i]*g.host[e][j]);
    std::array<Real,2>b={host[0][1],host[0][2]},center;
    for(int e=0;e<3;++e)for(int j=0;j<2;++j)center[j]=center[j]+Real(g.host[e][0]*g.host[e][1+j])/3;
    Real q=host[1][1]+host[2][2]-host[1][2]-gram_norm(b);static const int ij[3][2]={{0,1},{0,2},{1,2}};
    for(int active=1;active<=3;++active){auto&c=out[active-1];
        if(original){Real m[3][3];for(int i=0;i<2;++i)if(active&(1<<i))for(int e=0;e<3;++e)
            for(int j=0;j<3;++j)for(int k=0;k<3;++k)m[j][k]=m[j][k]+y[i]*child[i][ij[e][0]]*child[i][ij[e][1]]*Real(g.delta[i][e][j]*g.delta[i][e][k]);
            std::array<Real,2>cross={m[0][1],m[0][2]};Real energy=m[1][1]+m[2][2]-m[1][2];
            c.q=q;c.h=m[0][0];c.de=energy-2*gram_dot(b,cross)+c.h*gram_norm(b);c.ge=c.h*energy-gram_norm(cross);
            c.mean=energy-2*gram_dot(center,cross)+c.h*gram_norm(center)+c.h/12;
        }else{auto u=polynomial(g.U,alpha,child,y,active,false),f=polynomial(g.F,alpha,child,y,active,false),mean=polynomial(g.F,alpha,child,y,active,true);
            c.q=f[0];c.h=u[1];c.de=f[1]-c.h*c.q;c.ge=f[2];c.mean=mean[1]-c.h/4;}
    }return out;
}
Jet forest(const std::array<Coeff,3>&c,const Jet&mu,bool control){return primitive(c[2],mu,control)-primitive(c[0],mu,control)-primitive(c[1],mu,control);}
Real density(const Graph&g,bool original,bool control,const std::array<Real,8>&p){
    Real pi;arb_const_pi(pi.v,precision);Real r=root(1-p[0]);std::array<Real,3>alpha={p[0]/(1+r),r*p[1],r*(1-p[1])};
    std::array<std::array<Real,3>,2>child;Real jac(1);
    for(int i=0;i<2;++i){Real u=p[2+2*i],v=p[3+2*i],si,co,theta=pi*v/2;arb_sin_cos(si.v,co.v,theta.v,precision);
        Real rem=(1-u)*(1+u);child[i]={sq(u),rem*sq(si),rem*sq(co)};jac=jac*2*pi*u*rem*si*co;}
    Real si,co,theta=pi*p[6]/2;arb_sin_cos(si.v,co.v,theta.v,precision);std::array<Real,2>y={sq(si),sq(co)};
    auto data=coefficients(g,original,alpha,child,y);Real s=p[7],a=(1-s)*(1+s)/4,x=child[0][0];
    Real C=(a-2)/sq(1-a),D=1/(1-a),E=1/sq(1-a),mu=1+y[0]*x*(1/a-1);
    Jet unit=forest(data,Jet(Real(1),Real(1)),control),shifted=forest(data,Jet(mu),control);
    Real complete=C*unit.v+(E/a)*shifted.v-D*x*y[0]*unit.d;
    Real factor=pi*sq(sq(s))*jac/(24*si*co),answer=factor*complete;
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite complete spectral forest required");return answer;
}
Graph recipe(const char*path){std::ifstream input(path);std::string id,sha,extra;if(!(input>>id>>sha)||(input>>extra))throw std::invalid_argument("one spectral recipe required");auto g=load_graph(id);if(g.sha!=sha)throw std::invalid_argument("original graph digest required");return g;}
}
int main(int argc,char**argv){try{using namespace nested_spectral;if(argc!=2)throw std::invalid_argument("one recipe required");auto g=recipe(argv[1]);
    std::string mode;unsigned bits,control,count;if(!(std::cin>>mode>>bits>>control>>count)||(mode!="original"&&mode!="positive")||bits<512||bits>2048||control>1||count<1||count>512)throw std::invalid_argument("bounded spectral request required");
    precision=bits;std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"control\":"<<control<<",\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){std::array<Real,8>p;for(auto&v:p){auto q=read_rational();if(!(Rat(0)<q&&q<Rat(1)))throw std::invalid_argument("strict open eight-cube required");v=real(q);}if(i)std::cout<<',';interval(density(g,mode=="original",bool(control),p));}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing spectral request");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
