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
    std::vector<std::pair<int,int>> U,F,fullU,fullF; int delta[3][5]{};
    int gamma[64]{},nu[64]{},fraction[6]{};
    Rat weights[64][7],normalization;
    std::vector<Stratum> strata;
    int weighted_group=0;
};
#include "primitive_host_triangle_data.hpp"
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
// Full forest kernel. The unchanged prefix supplies exact Rat/CDF and Arb types.
Real triangle_log1p(const Real&x){Real y;arb_log1p(y.v,x.v,precision);return y;}
Real logarithm(const Real&x){Real y;arb_log(y.v,x.v,precision);return y;}
Real exponential(const Real&x){Real y;arb_exp(y.v,x.v,precision);return y;}
Real exprel(const Real&x){
    if(!arb_contains_zero(x.v)){Real y;arb_expm1(y.v,x.v,precision);return y/x;}
    // exp(convex hull of 0 and x) encloses integral_0^1 exp(t*x)dt, including x=0.
    Real y,zero(0);arb_union(y.v,x.v,zero.v,precision);arb_exp(y.v,y.v,precision);return y;
}
Real gramnorm(const std::array<Real,2>&v){return sq(v[0]-v[1]/2)+3*sq(v[1])/4;}
Gaussian matrix_gaussian(Real m[5][5]){
    Gaussian z;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){std::vector<int> r,c;for(int k=0;k<3;++k){if(k!=j)r.push_back(k);if(k!=i)c.push_back(k);}
        z.adj[i][j]=Real((i+j)%2?-1:1)*(m[r[0]][c[0]]*m[r[1]][c[1]]-m[r[0]][c[1]]*m[r[1]][c[0]]);}
    for(int i=0;i<3;++i)z.U=z.U+m[0][i]*z.adj[i][0];pos(z.U,"positive updated determinant required");
    Real reduced[2][2];for(int i=0;i<2;++i)for(int j=0;j<2;++j){reduced[i][j]=m[3+i][3+j];
        for(int r=0;r<3;++r)for(int c=0;c<3;++c)reduced[i][j]=reduced[i][j]-m[3+i][r]*z.adj[r][c]*m[c][3+j]/z.U;}
    z.Q=reduced[0][0]+reduced[1][1]-(reduced[0][1]+reduced[1][0])/2;return z;
}
Gaussian updated(const Graph&g,const std::array<Real,9>&alpha,const Real&r,const std::array<Real,3>&weights){
    Real m[5][5];for(int e=0;e<9;++e)for(int i=0;i<5;++i)for(int j=0;j<5;++j)m[i][j]=m[i][j]+alpha[e]*Real(long(g.routes[e][i]*g.routes[e][j]));
    for(int e=0;e<3;++e)for(int i=0;i<5;++i)for(int j=0;j<5;++j)m[i][j]=m[i][j]+r*weights[e]*Real(long(g.delta[e][i]*g.delta[e][j]));
    return matrix_gaussian(m);
}
Real gaussian_mean(const Graph&g,const std::array<Real,6>&a,const Gaussian&point,const Gaussian&center){
    Real q=center.Q;
    for(int j=0;j<6;++j)if(g.groups[j].size()==2){int e=g.groups[j][0];Real covariance;
        for(int k=0;k<3;++k)for(int l=0;l<3;++l)covariance=covariance+Real(long(g.routes[e][k]*g.routes[e][l]))*point.adj[k][l]/point.U;
        q=q-sq(a[j])*covariance/12;}
    return q;
}
template<size_t N>Real mean_polynomial(const Graph&g,const std::vector<std::pair<int,int>>&rows,const std::array<Real,6>&a,const std::array<Real,N>&alpha,int virtual_degree=-1){
    Real total;for(const auto&term:rows){int countvirtual=0;for(size_t e=9;e<N;++e)if(term.first&(1<<e))++countvirtual;
        if(virtual_degree>=0&&countvirtual!=virtual_degree)continue;Real monomial(term.second);
        for(int j=0;j<6;++j){int count=0;for(int e:g.groups[j])if(term.first&(1<<e))++count;
            if(g.groups[j].size()==2){if(count==1)monomial=monomial*a[j]/2;else if(count==2)monomial=monomial*sq(a[j])/6;}
            else if(count)monomial=monomial*a[j];}
        for(size_t e=9;e<N;++e)if(term.first&(1<<e))monomial=monomial*alpha[e];total=total+monomial;}
    return total;
}
Real linear_polynomial(const std::vector<std::pair<int,int>>&rows,const std::array<Real,12>&alpha){
    Real total;for(const auto&term:rows){int count=0;for(int e=9;e<12;++e)if(term.first&(1<<e))++count;if(count!=1)continue;
        Real monomial(term.second);for(int e=0;e<12;++e)if(term.first&(1<<e))monomial=monomial*alpha[e];total=total+monomial;}return total;
}
struct Boundary{Real ell,De,meanDe;};
Boundary boundary_original(const Graph&g,const Gaussian&host,const std::array<Real,9>&alpha,const std::array<Real,9>&central,const std::array<Real,6>&a,const std::array<Real,3>&w,bool control){
    Real B[3][2],Bc[3][2];for(int e=0;e<9;++e)for(int i=0;i<3;++i)for(int j=0;j<2;++j){
        B[i][j]=B[i][j]+alpha[e]*Real(long(g.routes[e][i]*g.routes[e][3+j]));Bc[i][j]=Bc[i][j]+central[e]*Real(long(g.routes[e][i]*g.routes[e][3+j]));}
    Boundary b;
    for(int d=0;d<3;++d){Real covariance;std::array<Real,2>res,cen;
        for(int j=0;j<2;++j){res[j]=cen[j]=Real(g.delta[d][3+j]);for(int k=0;k<3;++k)for(int l=0;l<3;++l){
            Real co=Real(g.delta[d][k])*host.adj[k][l]/host.U;res[j]=res[j]-co*B[l][j];cen[j]=cen[j]-co*Bc[l][j];}}
        for(int k=0;k<3;++k)for(int l=0;l<3;++l)covariance=covariance+Real(long(g.delta[d][k]*g.delta[d][l]))*host.adj[k][l]/host.U;
        b.ell=b.ell+w[d]*covariance;b.De=b.De+w[d]*gramnorm(res);
        if(control){b.meanDe=b.meanDe+w[d]*gramnorm(cen);
            for(int j=0;j<6;++j)if(g.groups[j].size()==2){int e=g.groups[j][0];Real z;
                for(int k=0;k<3;++k)for(int l=0;l<3;++l)z=z+Real(long(g.delta[d][k]*g.routes[e][l]))*host.adj[k][l]/host.U;
                b.meanDe=b.meanDe+w[d]*sq(a[j])*sq(z)/12;}}
    }
    return b;
}
Real triangle_density(const Graph&g,bool original,bool control,const Decoded&decoded,const std::array<Rat,11>&point){
    std::array<Real,6>a;Real cumulative(1),tropical(1);a[decoded.order[0]]=Real(1);
    for(int j=0;j<5;++j){Real ratio=root(real(point[j+1]),decoded.gamma[j]);cumulative=cumulative*ratio;a[decoded.order[j+1]]=cumulative;tropical=tropical*Pow::power(ratio,decoded.nu[j]);}
    std::array<Real,3>x={real(point[6]),real(point[7]),real(point[8])},xc={Real(1)/2,Real(1)/2,Real(1)/2};
    auto alpha=alphas(g,a,x),central=alphas(g,a,xc);Real M;for(const auto&v:a)M=M+v;
    auto host=gaussian(g,alpha);Real U0=original?host.U:poly(g.U,a),q0=original?host.Q:poly(g.F,alpha)/U0,mean0;
    if(control)mean0=original?gaussian_mean(g,a,host,gaussian(g,central)):mean_polynomial(g,g.F,a,alpha)/U0;
    Real first=real(point[9]),second=real(point[10]),s=root(1-first,2);std::array<Real,3>simplex={first/(1+s),s*second,s*(1-second)};
    static const int permutations[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    Real sum;
    for(const auto&permutation:permutations){std::array<Real,3>child={simplex[permutation[0]],simplex[permutation[1]],simplex[permutation[2]]};
        std::array<Real,3>w={child[0]*child[1],child[0]*child[2],child[1]*child[2]};
        std::array<Real,12>base;for(int e=0;e<9;++e)base[e]=alpha[e];for(int e=0;e<3;++e)base[9+e]=w[e];
        Boundary b;
        if(original)b=boundary_original(g,host,alpha,central,a,w,control);
        else{Real U1=linear_polynomial(g.fullU,base);b.ell=U1/U0;b.De=(linear_polynomial(g.fullF,base)-q0*U1)/U0;
            if(control)b.meanDe=(mean_polynomial(g,g.fullF,a,base,1)-mean0*U1)/U0;}
        pos(b.ell,"positive triangle loop variance required");Real h=logarithm(M*b.ell);
        for(int branch=0;branch<2;++branch){Rat uq=branch?Rat(1)-decoded.auxiliary:decoded.auxiliary;Real u=real(uq),contribution;
            if(fmpq_is_zero(uq.v)){
                contribution=M/exprel(h)*(3*b.ell*triangle_log1p(q0/M)-b.De/(M+q0));
                if(control)contribution=contribution+M/exprel(h)*(3*b.ell*(mean0-q0)-(b.meanDe-b.De))/M;
            }else if(fmpq_is_one(uq.v))contribution=(control?mean0:q0)/(M*exprel(-h));
            else{
                Real eu=exprel(u*h),ev=exprel((1-u)*h);
                Real r=M*u/(1-u)*eu/(exponential(u*h)*ev),jac=exprel(h)/(u*(1-u)*eu*ev),mass=M+r;
                Real Ur,qr,meanr;
                if(original){auto current=updated(g,alpha,r,w);Ur=current.U;qr=current.Q;
                    if(control)meanr=gaussian_mean(g,a,current,updated(g,central,r,w));}
                else{auto active=base;for(int e=9;e<12;++e)active[e]=r*active[e];Ur=poly(g.fullU,active);qr=poly(g.fullF,active)/Ur;
                    if(control)meanr=mean_polynomial(g,g.fullF,a,active)/Ur;}
                Real determinant_delta=Ur/U0-1,ratio=1+determinant_delta,inv3=1/Pow::power(ratio,3),delta=qr-q0;
                Real one_minus=determinant_delta*(3+3*determinant_delta+sq(determinant_delta))*inv3;
                contribution=one_minus*triangle_log1p(q0/mass)-inv3*triangle_log1p(delta/(mass+q0));
                if(control)contribution=contribution+(one_minus*(mean0-q0)-inv3*((meanr-mean0)-delta))/mass;
                contribution=contribution*jac;
            }
            sum=sum+contribution;
        }
    }
    return real(g.normalization)*Pow::power(tropical/U0,3)*sum/24;
}
Real triangle_value(const Graph&g,bool original,bool control,const std::array<Rat,11>&point){
    Real total;for(const auto&st:g.strata)total=total+real(st.probability)*triangle_density(g,original,control,decode(g,st,point[0]),point);
    if(!arb_is_finite(total.v))throw std::runtime_error("finite complete triangle density required");return total;
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
        throw std::invalid_argument("bounded triangle host request required");
    precision=bits;
    std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"control\":"<<control
             <<",\"strata\":26,\"radial_reflection\":2,\"simplex_permutations\":6,\"matched_radial_transport\":1,\"count\":"<<count<<",\""<<(mode=="decode"?"decodings":"values")<<"\":[";
    for(unsigned i=0;i<count;++i){std::array<Rat,11>p;for(auto&q:p){q=read_rational();if(!(Rat(0)<q)||!(q<Rat(1)))throw std::invalid_argument("strict open eleven-cube required");}
        if(i)std::cout<<',';
        if(mode=="decode"){std::cout<<'[';for(size_t j=0;j<g.strata.size();++j){if(j)std::cout<<',';print_decoded(decode(g,g.strata[j],p[0]));}std::cout<<']';}
        else interval(triangle_value(g,mode=="original",bool(control),p));}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing triangle request data");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
