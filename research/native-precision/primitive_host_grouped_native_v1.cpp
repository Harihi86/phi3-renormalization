// Seven complete four-loop graphs: exact discrete decoder and Arb integrands.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include "fmpq.h"
namespace grouped_host {
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::interval;
namespace Pow=disjoint_v_integrated;
struct Rat {
    fmpq_t v;
    Rat(){fmpq_init(v);}
    explicit Rat(long x){fmpq_init(v);fmpq_set_si(v,x,1);}
    explicit Rat(const std::string&s){fmpq_init(v);if(fmpq_set_str(v,s.c_str(),10)||fmpz_sgn(fmpq_denref(v))<=0){fmpq_clear(v);throw std::invalid_argument("rational with positive denominator required");}fmpq_canonicalise(v);}
    Rat(const Rat&o){fmpq_init(v);fmpq_set(v,o.v);}
    Rat& operator=(const Rat&o){fmpq_set(v,o.v);return *this;}
    ~Rat(){fmpq_clear(v);}
    std::string str()const{char*p=fmpq_get_str(nullptr,10,v);std::string s(p);flint_free(p);return s;}
};
Rat operator+(const Rat&a,const Rat&b){Rat z;fmpq_add(z.v,a.v,b.v);return z;}
Rat operator-(const Rat&a,const Rat&b){Rat z;fmpq_sub(z.v,a.v,b.v);return z;}
Rat operator*(const Rat&a,const Rat&b){Rat z;fmpq_mul(z.v,a.v,b.v);return z;}
Rat operator/(const Rat&a,const Rat&b){if(fmpq_is_zero(b.v))throw std::runtime_error("zero rational divisor");Rat z;fmpq_div(z.v,a.v,b.v);return z;}
bool operator<(const Rat&a,const Rat&b){return fmpq_cmp(a.v,b.v)<0;}
Real real(const Rat&q){Real z;arb_set_fmpq(z.v,q.v,precision);return z;}
Rat read_rational(){std::string n,d;if(!(std::cin>>n>>d)||n.size()>10000||d.size()>10000)throw std::invalid_argument("bounded rational input required");return Rat(n+"/"+d);}
struct Stratum{int first,position;Rat probability;};
struct Graph {
    std::string id,sha;
    int routes[9][5]{};
    std::vector<std::vector<int>> groups;
    std::vector<std::pair<int,int>> U,F;
    int gamma[64]{},nu[64]{},fraction[6]{};
    Rat weights[64][7],normalization;
    std::vector<Stratum> strata;
    int weighted_group=0;
};
#include "primitive_host_grouped_data.hpp"
struct Decoded{std::array<int,6> order{};std::array<int,5> gamma{},nu{};Rat auxiliary,probability{1};};
Decoded decode(const Graph&g,const Stratum&s,Rat z){
    if(z<Rat(0)||!(z<Rat(1)))throw std::invalid_argument("unit-interval CDF point required");
    Decoded out;int remaining=63;
    for(int depth=0;depth<6;++depth){
        int k=(s.position!=-1&&(remaining&(1<<g.weighted_group)))?s.position-depth:-1;
        std::vector<int> options,children;std::vector<Rat> masses;Rat total;
        for(int e=0;e<6;++e)if((remaining&(1<<e))&&(k==-1||((e==g.weighted_group)==(k==0)))&&(depth!=0||s.first==-1||e==s.first)){
            int child=remaining^(1<<e),next=(k==-1||e==g.weighted_group)?-1:k-1;
            Rat mass=g.weights[child][next+1]/Rat(child?g.gamma[child]:1);
            options.push_back(e);children.push_back(child);masses.push_back(mass);total=total+mass;
        }
        Rat left;bool selected=false;
        for(size_t j=0;j<options.size();++j){Rat p=masses[j]/total;
            if(z<left+p){z=(z-left)/p;out.probability=out.probability*p;out.order[depth]=options[j];remaining=children[j];selected=true;break;}
            left=left+p;
        }
        if(!selected)throw std::runtime_error("exact conditional CDF failed");
        if(depth<5){out.gamma[depth]=g.gamma[remaining];out.nu[depth]=g.nu[remaining];}
    }
    if(remaining||out.order[0]!=s.first||out.order[s.position]!=g.weighted_group)throw std::runtime_error("sector constraint failed");
    out.auxiliary=z;return out;
}
Real sq(const Real&x){Real z;arb_sqr(z.v,x.v,precision);return z;}
Real root(const Real&x,unsigned n){Real z;arb_root_ui(z.v,x.v,n,precision);return z;}
void pos(const Real&x,const char*s){if(!arb_is_positive(x.v))throw std::runtime_error(s);}
Real bubble(const Real&t,bool original){
    pos(t,"positive bubble argument required");
    if(original){auto z=disjoint_interval::bubble(fourier_e2_root::Z(t));if(!arb_contains_zero(acb_imagref(z.v)))throw std::runtime_error("bubble is not real");Real r;arb_set(r.v,acb_realref(z.v));return r;}
    Real r=root(t/(t+4),2),v;pos(1-r,"atanh interior required");arb_atanh(v.v,r.v,precision);
    return (t+4)/3*(v/r-1)-t/9;
}
template<size_t N>Real poly(const std::vector<std::pair<int,int>>&rows,const std::array<Real,N>&a){
    Real sum;for(const auto&entry:rows){Real monomial(entry.second);for(size_t e=0;e<N;++e)if(entry.first&(1<<e))monomial=monomial*a[e];sum=sum+monomial;}return sum;
}
std::array<Real,9> alphas(const Graph&g,const std::array<Real,6>&a,const std::array<Real,3>&x){
    std::array<Real,9> result;
    for(int j=0;j<6;++j){const auto&edges=g.groups[j];if(edges.size()==1)result[edges[0]]=a[j];
        else{const auto&v=x[g.fraction[j]];result[edges[0]]=a[j]*v;result[edges[1]]=a[j]*(1-v);}}
    return result;
}
struct Gaussian{Real U,Q;Real adj[3][3];};
Gaussian gaussian(const Graph&g,const std::array<Real,9>&alpha){
    Real m[5][5];for(int e=0;e<9;++e)for(int i=0;i<5;++i)for(int j=0;j<5;++j)
        m[i][j]=m[i][j]+alpha[e]*Real(long(g.routes[e][i]*g.routes[e][j]));
    Gaussian z;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){std::vector<int> r,c;for(int k=0;k<3;++k){if(k!=j)r.push_back(k);if(k!=i)c.push_back(k);}
        z.adj[i][j]=Real((i+j)%2?-1:1)*(m[r[0]][c[0]]*m[r[1]][c[1]]-m[r[0]][c[1]]*m[r[1]][c[0]]);}
    for(int i=0;i<3;++i)z.U=z.U+m[0][i]*z.adj[i][0];pos(z.U,"positive original determinant required");
    Real reduced[2][2];for(int i=0;i<2;++i)for(int j=0;j<2;++j){reduced[i][j]=m[3+i][3+j];
        for(int r=0;r<3;++r)for(int c=0;c<3;++c)reduced[i][j]=reduced[i][j]-m[3+i][r]*z.adj[r][c]*m[c][3+j]/z.U;}
    z.Q=reduced[0][0]+reduced[1][1]-(reduced[0][1]+reduced[1][0])/2;return z;
}
Real density(const Graph&g,bool original,bool control,const Decoded&decoded,const std::array<Rat,9>&point){
    std::array<Real,6>a;Real cumulative(1),tropical(1);a[decoded.order[0]]=Real(1);
    for(int j=0;j<5;++j){Real ratio=root(real(point[j+1]),decoded.gamma[j]);cumulative=cumulative*ratio;
        a[decoded.order[j+1]]=cumulative;tropical=tropical*Pow::power(ratio,decoded.nu[j]);}
    std::array<Real,3>x={real(point[6]),real(point[7]),real(point[8])};Real jac(1),M;
    for(int j=0;j<6;++j){M=M+a[j];if(g.groups[j].size()==2&&j==g.weighted_group){int k=g.fraction[j];jac=2*x[k];x[k]=sq(x[k]);}}
    auto alpha=alphas(g,a,x);Real U,Q,Qmean;
    if(original){auto exact=gaussian(g,alpha);U=exact.U;Q=exact.Q;
        if(control){auto center=x;for(int j=0;j<6;++j)if(g.groups[j].size()==2&&j!=g.weighted_group)center[g.fraction[j]]=Real(1)/2;
            Qmean=gaussian(g,alphas(g,a,center)).Q;
            for(int j=0;j<6;++j)if(g.groups[j].size()==2&&j!=g.weighted_group){int edge=g.groups[j][0];Real correction;
                for(int r=0;r<3;++r)for(int c=0;c<3;++c)correction=correction+Real(long(g.routes[edge][r]*g.routes[edge][c]))*exact.adj[r][c];
                Qmean=Qmean-sq(a[j])*correction/(12*U);}}
    }else{
        U=poly(g.U,a);Q=poly(g.F,alpha)/U;
        if(control){Real meanF;
            for(const auto&term:g.F){Real monomial(term.second);
                for(int j=0;j<6;++j){const auto&edges=g.groups[j];int count=0;for(int e:edges)if(term.first&(1<<e))++count;
                    if(edges.size()==2&&j!=g.weighted_group){if(count==1)monomial=monomial*a[j]/2;else if(count==2)monomial=monomial*sq(a[j])/6;}
                    else for(int e:edges)if(term.first&(1<<e))monomial=monomial*alpha[e];}
                meanF=meanF+monomial;}Qmean=meanF/U;}
    }
    pos(U,"positive U required");pos(Q,"positive physical external form required");
    Real d=alpha[8],other;for(int e=0;e<8;++e)other=other+alpha[e];
    Real u=real(decoded.auxiliary),value;
    for(int branch=0;branch<2;++branch){Real v=branch?1-u:u,mass=M+v*Q;
        value=value-d*Q*bubble((other+v*Q)/d,original)/(2*sq(mass));}
    if(control)value=value+d*bubble(other/d,original)/sq(M)*(Q-Qmean);
    return real(g.normalization)*Pow::power(tropical/U,3)*jac*value;
}
Real value(const Graph&g,bool original,bool control,const std::array<Rat,9>&p){
    Real total;for(const auto&s:g.strata)total=total+real(s.probability)*density(g,original,control,decode(g,s,p[0]),p);
    if(!arb_is_finite(total.v))throw std::runtime_error("finite complete density required");return total;
}
void print_decoded(const Decoded&d){
    std::cout<<"{\"order\":[";for(int j=0;j<6;++j){if(j)std::cout<<',';std::cout<<d.order[j];}
    std::cout<<"],\"exponents\":[";for(int j=0;j<5;++j){if(j)std::cout<<',';std::cout<<d.gamma[j];}
    std::cout<<"],\"leading\":[";for(int j=0;j<5;++j){if(j)std::cout<<',';std::cout<<d.nu[j];}
    std::cout<<"],\"auxiliary\":\""<<d.auxiliary.str()<<"\",\"probability\":\""<<d.probability.str()<<"\"}";
}
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
    if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("exact graph recipe required");Graph g=load_graph(id);
    if(sha!=g.sha||g.strata.size()!=26)throw std::invalid_argument("frozen original graph required");return g;}
}
int main(int argc,char**argv){try{using namespace grouped_host;
    if(argc!=2)throw std::invalid_argument("one recipe path required");Graph g=recipe(argv[1]);
    std::string mode;unsigned bits,control,count;
    if(!(std::cin>>mode>>bits>>control>>count)||(mode!="closed"&&mode!="original"&&mode!="decode")||bits<512||bits>2048||control>1||count<1||count>512)
        throw std::invalid_argument("bounded grouped-host request required");
    precision=bits;
    std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"control\":"<<control
             <<",\"strata\":26,\"auxiliary_antithetic\":2,\"weighted_fraction_power\":2,\"count\":"<<count<<",\""<<(mode=="decode"?"decodings":"values")<<"\":[";
    for(unsigned i=0;i<count;++i){std::array<Rat,9>p;for(auto&q:p){q=read_rational();if(!(Rat(0)<q)||!(q<Rat(1)))throw std::invalid_argument("strict open nine-cube required");}
        if(i)std::cout<<',';
        if(mode=="decode"){std::cout<<'[';for(size_t j=0;j<g.strata.size();++j){if(j)std::cout<<',';print_decoded(decode(g,g.strata[j],p[0]));}std::cout<<']';}
        else interval(value(g,mode=="original",bool(control),p));}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing request data");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
