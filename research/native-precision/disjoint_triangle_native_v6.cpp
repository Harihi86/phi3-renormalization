// Complete disjoint-triangle forests on a common radial chart.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include "fmpq.h"
namespace disjoint_triangle {
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
struct Graph{std::string id,sha;int L=0,k=0,count=0;std::vector<std::vector<int>> routes,groups;std::vector<std::vector<std::vector<int>>> delta;std::vector<std::pair<int,int>> U,F;};
#include "disjoint_triangle_data.hpp"
using Cube=std::vector<Real>;
using Matrix=std::vector<std::vector<Cube>>;
Real sq(const Real&x){Real y;arb_sqr(y.v,x.v,precision);return y;}
Real root(const Real&x){Real y;arb_sqrt(y.v,x.v,precision);return y;}
void pos(const Real&x,const char*s){if(!arb_is_positive(x.v))throw std::runtime_error(s);}
Real logarithm(const Real&x){Real y;arb_log(y.v,x.v,precision);return y;}
Real stable_log1p(const Real&x){Real y;arb_log1p(y.v,x.v,precision);return y;}
Real exponential(const Real&x){Real y;arb_exp(y.v,x.v,precision);return y;}
Real exprel(const Real&x){
    if(!arb_contains_zero(x.v)){Real y;arb_expm1(y.v,x.v,precision);return y/x;}
    Real y,zero(0);arb_union(y.v,x.v,zero.v,precision);arb_exp(y.v,y.v,precision);return y;
}
Real divided_log(const Real&t,const Real&positive,const Real&negative){
    if(arb_contains_zero(t.v)){Real h,zero(0);arb_union(h.v,t.v,zero.v,precision);pos(1+h,"positive divided-log hull");return 1/(1+h);}
    Real ab;arb_abs(ab.v,t.v);Real quarter=Real(1)/4;
    if(arb_lt(ab.v,quarter.v))return stable_log1p(t)/t;
    pos(positive,"positive logarithm product");pos(negative,"negative-parity product is positive");
    return (logarithm(positive)-logarithm(negative))/t;
}
int pop(int mask){return __builtin_popcount(unsigned(mask));}
Cube constant(const Real&value,int k){Cube out(1<<k);out[0]=value;return out;}
Cube add(const Cube&a,const Cube&b){Cube out(a.size());for(size_t i=0;i<a.size();++i)out[i]=a[i]+b[i];return out;}
Cube sub(const Cube&a,const Cube&b){Cube out(a.size());for(size_t i=0;i<a.size();++i)out[i]=a[i]-b[i];return out;}
Cube scale(const Cube&a,const Real&v){Cube out(a.size());for(size_t i=0;i<a.size();++i)out[i]=a[i]*v;return out;}
struct Algebra{
    int k,n;std::vector<Real> r,w;
    explicit Algebra(const std::vector<Real>&v):k(int(v.size())),n(1<<k),r(v),w(n){
        w[0]=Real(1);for(int m=1;m<n;++m){int i=__builtin_ctz(unsigned(m));w[m]=w[m^(1<<i)]*r[i];}}
    Cube mul(const Cube&a,const Cube&b)const{Cube out(n);for(int i=0;i<n;++i)if(!arb_is_zero(a[i].v))for(int j=0;j<n;++j)if(!arb_is_zero(b[j].v))out[i|j]=out[i|j]+a[i]*b[j]*w[i&j];return out;}
    Real vertex(const Cube&a,int mask)const{Real out;for(int i=0;i<n;++i)if((i&mask)==i)out=out+a[i]*w[i];return out;}
    Cube restrict(const Cube&a,int mask)const{Cube out=a;for(int i=0;i<n;++i)if((i&mask)!=i)out[i]=Real(0);return out;}
    Cube inverse(const Cube&a)const{Cube out(n);pos(a[0],"positive host cube inverse");out[0]=1/a[0];
        for(int degree=1;degree<=k;++degree)for(int mask=1;mask<n;++mask)if(pop(mask)==degree){Real numerator;
            for(int b=0;b<n;++b)if(b!=mask&&(b&mask)==b)for(int c=0;c<n;++c)if((c&mask)==c&&(c|b)==mask)numerator=numerator+a[c]*out[b]*w[c&b];
            Real denominator=vertex(a,mask);pos(denominator,"positive inverse cube corner");out[mask]=-numerator/denominator;}
        return out;}
    Cube power(const Cube&a,int degree)const{Cube out=constant(Real(1),k);for(int i=0;i<degree;++i)out=mul(out,a);return out;}
    Cube log(const Cube&a)const{Cube out(n);pos(a[0],"positive cube logarithm");out[0]=logarithm(a[0]);
        for(int mask=1;mask<n;++mask){Cube prod[2]={constant(Real(1),k),constant(Real(1),k)};Real value[2]={Real(1),Real(1)};
            for(int part=0;part<n;++part)if((part&mask)==part){int parity=(pop(mask)-pop(part))%2;prod[parity]=mul(prod[parity],restrict(a,part));value[parity]=value[parity]*vertex(a,part);}
            pos(value[1],"positive cube log denominator");Real factor=(prod[0][mask]-prod[1][mask])/value[1],t=w[mask]*factor;
            out[mask]=factor*divided_log(t,value[0],value[1]);}
        return out;}
};
Cube determinant(const Matrix&m,const Algebra&A){int n=int(m.size());std::vector<int>p(n);for(int i=0;i<n;++i)p[i]=i;Cube out=constant(Real(0),A.k);
    do{int inversions=0;Cube term=constant(Real(1),A.k);for(int i=0;i<n;++i){for(int j=i+1;j<n;++j)if(p[i]>p[j])++inversions;term=A.mul(term,m[i][p[i]]);}out=add(out,scale(term,Real(inversions%2?-1:1)));}while(std::next_permutation(p.begin(),p.end()));return out;}
Matrix minor(const Matrix&m,const std::vector<int>&rows,const std::vector<int>&cols){Matrix out;for(int i:rows){std::vector<Cube>r;for(int j:cols)r.push_back(m[i][j]);out.push_back(r);}return out;}
Matrix gaussian(const Graph&g,const std::vector<Real>&alpha,const std::vector<std::array<Real,3>>&child,const Algebra&A){
    Matrix out(g.L+2,std::vector<Cube>(g.L+2,constant(Real(0),g.k)));static const int ij[3][2]={{0,1},{0,2},{1,2}};
    for(int a=0;a<g.L+2;++a)for(int b=0;b<g.L+2;++b){
        for(int e=0;e<g.count;++e)out[a][b][0]=out[a][b][0]+alpha[e]*Real(g.routes[e][a]*g.routes[e][b]);
        for(int c=0;c<g.k;++c)for(int j=0;j<3;++j)out[a][b][1<<c]=out[a][b][1<<c]+child[c][ij[j][0]]*child[c][ij[j][1]]*Real(g.delta[c][j][a]*g.delta[c][j][b]);}
    return out;}
std::pair<Cube,Cube> gaussian_polynomials(const Graph&g,const Matrix&m,const Algebra&A){std::vector<int>loop;for(int j=0;j<g.L;++j)loop.push_back(j);
    Cube U=determinant(minor(m,loop,loop),A),F=constant(Real(0),g.k);
    for(int i=0;i<2;++i)for(int j=0;j<2;++j){auto rows=loop,cols=loop;rows.push_back(g.L+i);cols.push_back(g.L+j);F=add(F,scale(determinant(minor(m,rows,cols),A),i==j?Real(1):-Real(1)/2));}
    return {U,F};}
Cube polynomial(const Graph&g,const std::vector<std::pair<int,int>>&rows,const std::vector<Real>&alpha,const std::vector<std::array<Real,3>>&child,const Algebra&A,const std::array<Real,3>&a,bool mean){
    static const int ij[3][2]={{0,1},{0,2},{1,2}};Cube out=constant(Real(0),g.k);
    for(const auto&entry:rows){int mask=entry.first,bits=0;Real term(entry.second);
        if(!mean){for(int e=0;e<g.count;++e)if(mask&(1<<e))term=term*alpha[e];}
        else if(g.L==1){int degree=pop(mask&7),factorial=1;for(int j=2;j<=degree+2;++j)factorial*=j;term=term*2/factorial;}
        else for(int j=0;j<3;++j){int degree=0;for(int e:g.groups[j])if(mask&(1<<e))++degree;if(degree==1)term=term*a[j]/2;else if(degree==2)term=term*sq(a[j])/6;}
        for(int c=0;c<g.k;++c){int degree=0;for(int j=0;j<3;++j)if(mask&(1<<(g.count+3*c+j))){++degree;term=term*child[c][ij[j][0]]*child[c][ij[j][1]];}
            if(degree){bits|=1<<c;for(int j=1;j<degree;++j)term=term*A.r[c];}}
        out[bits]=out[bits]+term;}
    return out;}
Cube gaussian_mean(const Graph&g,const Matrix&point,const Matrix&center,const Algebra&A,const std::array<Real,3>&a){
    Cube mean=gaussian_polynomials(g,center,A).second;
    if(g.L==1)mean[0]=mean[0]-Real(1)/12;
    else for(int j=0;j<3;++j){int e=g.groups[j][0],v=g.routes[e][0],w=g.routes[e][1];
        Cube form=add(scale(point[1][1],Real(v*v)),scale(point[0][0],Real(w*w)));form=sub(form,scale(add(point[0][1],point[1][0]),Real(v*w)));
        mean=sub(mean,scale(form,sq(a[j])/12));}
    return mean;}
Real forest(Cube U,Cube F,Cube mean,const Real&mass,const Algebra&A,bool control){Cube inverse=A.inverse(U);
    Cube logs=sub(A.log(add(U,scale(F,1/mass))),A.log(U));logs[0]=stable_log1p(F[0]/(mass*U[0]));
    Cube value=A.mul(A.power(inverse,3),logs);
    if(control)value=add(value,scale(A.mul(sub(mean,F),A.power(inverse,4)),1/mass));
    return -value.back();}
std::array<Real,3> simplex(const Real&u,const Real&v){Real s=root(1-u);return {u/(1+s),s*v,s*(1-v)};}
Real density(const Graph&g,bool original,bool control,int reflection,int transport,const std::array<Rat,11>&input){
    std::array<Real,11>p;for(int j=0;j<11;++j)p[j]=real(input[j]);std::array<int,3>order={0,1,2};Real total;int orders=g.L==2?6:1;
    do{std::array<Real,3>a;std::vector<Real>alpha(g.count),central(g.count);Real tropical(1),M;int base;
        if(g.L==2){a[order[0]]=Real(1);a[order[1]]=p[0];a[order[2]]=p[0]*root(p[1]);tropical=a[order[0]]*a[order[1]];base=5;
            for(int j=0;j<3;++j){int e=g.groups[j][0],f=g.groups[j][1];alpha[e]=a[j]*p[2+j];alpha[f]=a[j]*(1-p[2+j]);central[e]=central[f]=a[j]/2;}}
        else{auto x=simplex(p[0],p[1]);for(int e=0;e<3;++e){alpha[e]=x[e];central[e]=Real(1)/3;}base=2;}
        for(const auto&v:alpha)M=M+v;
        std::vector<std::array<Real,3>>child;for(int i=0;i<g.k;++i)child.push_back(simplex(p[base+2*i],p[base+2*i+1]));
        Real host[2][2];for(int e=0;e<g.count;++e)for(int i=0;i<g.L;++i)for(int j=0;j<g.L;++j)host[i][j]=host[i][j]+alpha[e]*Real(g.routes[e][i]*g.routes[e][j]);
        Real U0=g.L==1?host[0][0]:host[0][0]*host[1][1]-host[0][1]*host[1][0];pos(U0,"positive host determinant");
        std::vector<Real>shape;if(g.k==2)shape={p[9],1-p[9]};else{auto y=simplex(p[8],p[9]);shape={y[0],y[1],y[2]};}
        Real ell;std::vector<Real> childell(g.k);static const int ij[3][2]={{0,1},{0,2},{1,2}};
        for(int c=0;c<g.k;++c)for(int j=0;j<3;++j){int v=g.delta[c][j][0],w=g.L==2?g.delta[c][j][1]:0;Real cov;
            if(g.L==1)cov=Real(v*v)/U0;else cov=(Real(v*v)*host[1][1]+Real(w*w)*host[0][0]-Real(v*w)*(host[0][1]+host[1][0]))/U0;
            childell[c]=childell[c]+child[c][ij[j][0]]*child[c][ij[j][1]]*cov;}
        if(transport==2){Real denominator;for(int c=0;c<g.k;++c){pos(childell[c],"positive individual child variance");denominator=denominator+1/(shape[c]*childell[c]);}ell=1/denominator;}
        else for(int c=0;c<g.k;++c)ell=ell+shape[c]*childell[c];
        pos(ell,"positive common radial loop variance");Real h=logarithm(M*ell);
        for(int branch=0;branch<reflection;++branch){Real u=branch?1-p[10]:p[10],radial=M*u/(1-u),derivative=M/sq(1-u);
            if(transport){Real eu=exprel(u*h),ev=exprel((1-u)*h);radial=M*u/(1-u)*eu/(exponential(u*h)*ev);derivative=radial*exprel(h)/(u*(1-u)*eu*ev);}
            std::vector<Real>r;for(const auto&y:shape)r.push_back(radial*y);Algebra A(r);Cube U,F,mean;
            if(original){auto m=gaussian(g,alpha,child,A),centralmatrix=gaussian(g,central,child,A);auto polys=gaussian_polynomials(g,m,A);U=polys.first;F=polys.second;if(control)mean=gaussian_mean(g,m,centralmatrix,A,a);}
            else{U=polynomial(g,g.U,alpha,child,A,a,false);F=polynomial(g,g.F,alpha,child,A,a,false);if(control)mean=polynomial(g,g.F,alpha,child,A,a,true);}
            U=scale(U,1/U0);F=scale(F,1/U0);if(control)mean=scale(mean,1/U0);
            Real jac=Pow::power(radial,g.k-1)*derivative/(g.k==3?2:1),volume=Real(1)/((1<<g.k)*(g.L==1?2:1));
            total=total+Real(g.L==2?3:1)*Pow::power(tropical/U0,3)*forest(U,F,mean,M+radial,A,control)*jac*volume/(orders*reflection);
        }
        if(g.L==1)break;
    }while(std::next_permutation(order.begin(),order.end()));
    if(!arb_is_finite(total.v))throw std::runtime_error("finite common radial density required");return total;}
unsigned radial_series_terms(const Real&upper){
    if(arb_is_zero(upper.v))return 24;
    slong exponent=arf_abs_bound_lt_2exp_si(arb_midref(upper.v)),decay=std::max(slong(1),-exponent);
    return std::max(24u,unsigned((precision+decay-1)/decay)+8);
}
Real radial_phi(const Real&t){
    pos(t,"positive radial phi argument");Real pi;arb_const_pi(pi.v,precision);Real constant=sq(pi)/6-Real(5)/4;
    if(arb_is_one(t.v))return constant;
    Real delta=t-1,upper=arb_cw::abs_upper(delta),half=Real(1)/2;
    if(arb_le(upper.v,half.v)){
        unsigned count=radial_series_terms(upper);Real result=constant,power=delta;
        for(unsigned j=1;j<=count;++j){result=result+Real(j%2?2:-2)*power/(Real(long(j))*Real(long(j))*Real(long(j+1))*Real(long(j+2)));power=power*delta;}
        Real tail=2*Pow::power(upper,count+1)/(Real(long(count+1))*Real(long(count+1))*Real(long(count+2))*Real(long(count+3))*(1-upper));
        arb_cw::add_error(result,arb_cw::abs_upper(tail));return result;
    }
    Real dilog,arg=1-t;arb_polylog_si(dilog.v,2,arg.v,precision);
    return sq(pi)/6-dilog+t*(2-3*t)*logarithm(t)/(2*sq(1-t))-t/(2*(1-t));
}
std::array<Real,4> radial_moments(const Real&h){
    Real delta=h-1,upper=arb_cw::abs_upper(delta),half=Real(1)/2;std::array<Real,4>out;
    if(arb_le(upper.v,half.v)){
        unsigned count=radial_series_terms(upper);
        for(unsigned n=1;n<=4;++n){Real power(1),sum;for(unsigned j=0;j<=count;++j){sum=sum+power/Real(long(n+j));power=-power*delta;}
            Real tail=Pow::power(upper,count+1)/(Real(long(n+count+1))*(1-upper));arb_cw::add_error(sum,arb_cw::abs_upper(tail));out[n-1]=sum;}
    }else{out[0]=logarithm(h)/delta;for(int n=2;n<=4;++n)out[n-1]=(Real(1)/Real(n-1)-out[n-2])/delta;}
    return out;
}
Real radial_primitive(const Real&q,const Real&h,const Real&de,const Real&ge,const Real&mean_de,bool control){
    pos(h,"positive nonempty child variance");Real a=1+q,b=1+h*a+de,c=h+ge,variance=h*de-ge;
    // The exact Schur identity proves variance=|c_vector-h*b_vector|^2.
    // Intersecting its interval with the nonnegative axis is therefore valid.
    if(arb_is_negative(variance.v))throw std::runtime_error("inconsistent exact Schur variance");
    arb_nonnegative_part(variance.v,variance.v);
    Real disc=sq(1-h*a+de)+4*a*variance,s=root(disc),plus=(b+s)/(2*a),minus=2*c/(b+s),la=stable_log1p(q);
    Real result=la*logarithm(h)+sq(la)/2+3*la/2-radial_phi(plus/h)-radial_phi(minus/h)+radial_phi(1/h)+radial_phi(Real(1));
    if(control){auto moments=radial_moments(h);result=result+(Real(1)/4-q)*h*(moments[0]+moments[1]+moments[2])-(mean_de-de)*moments[3];}
    return result;
}
std::array<Real,3> common_radial_coefficients(const Graph&g,const std::vector<std::pair<int,int>>&rows,const std::vector<Real>&alpha,const std::vector<std::array<Real,3>>&child,const std::array<Real,3>&y,int active,bool mean){
    std::array<Real,3>out;static const int ij[3][2]={{0,1},{0,2},{1,2}};
    for(const auto&entry:rows){int mask=entry.first,degree=0;Real monomial(entry.second);bool use=true;
        if(mean){int d=pop(mask&7),factorial=1;for(int j=2;j<=d+2;++j)factorial*=j;monomial=monomial*2/factorial;}
        else for(int e=0;e<3;++e)if(mask&(1<<e))monomial=monomial*alpha[e];
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(mask&(1<<(3+3*i+j))){
            if(!(active&(1<<i))){use=false;break;}
            ++degree;monomial=monomial*y[i]*child[i][ij[j][0]]*child[i][ij[j][1]];}
        if(use){if(degree>2)throw std::runtime_error("unexpected radial polynomial degree");out[degree]=out[degree]+monomial;}
    }
    return out;
}
Real gram_norm(const std::array<Real,2>&v){return sq(v[0]-v[1]/2)+3*sq(v[1])/4;}
Real gram_dot(const std::array<Real,2>&v,const std::array<Real,2>&w){return (v[0]-v[1]/2)*(w[0]-w[1]/2)+3*v[1]*w[1]/4;}
Real radial_density(const Graph&g,bool original,bool control,const std::array<Rat,11>&input,int shape_mode=0,int permutation=0,int mixture_branch=-1){
    if(g.id!="L4E3G521"||g.L!=1||g.k!=3)throw std::invalid_argument("integrated radial mode requires the proved three-child graph");
    std::array<Real,11>p;for(int j=0;j<11;++j)p[j]=real(input[j]);auto hostx=simplex(p[0],p[1]),y=simplex(p[8],p[9]);
    std::vector<Real>alpha={hostx[0],hostx[1],hostx[2]};std::vector<std::array<Real,3>>child;
    Real childjac(1);
    for(int i=0;i<3;++i){
        if(shape_mode>=3){Real pi;arb_const_pi(pi.v,precision);Real u=p[2+2*i],v=p[3+2*i],theta=pi*v/2,sine,cosine;
            arb_sin_cos(sine.v,cosine.v,theta.v,precision);Real rem=(1-u)*(1+u);child.push_back({sq(u),rem*sq(sine),rem*sq(cosine)});
            childjac=childjac*2*pi*u*rem*sine*cosine;
        }else child.push_back(simplex(p[2+2*i],p[3+2*i]));}

    Real density_factor=Real(1)/(32*y[0]*y[1]*y[2]);
    if(shape_mode){
        Real pi;arb_const_pi(pi.v,precision);Real theta=pi*p[9]/2,sine,cosine;
        arb_sin_cos(sine.v,cosine.v,theta.v,precision);Real rem=(1-p[8])*(1+p[8]);
        y={sq(p[8]),rem*sq(sine),rem*sq(cosine)};
        // sqrt(prod w) is evaluated with positive factors, without subtraction.
        density_factor=pi/(8*p[8]*rem*sine*cosine);
        if(shape_mode>=3)density_factor=density_factor*8*childjac;
        if(permutation){std::array<int,3>order={0,1,2};for(int j=0;j<permutation;++j)std::next_permutation(order.begin(),order.end());auto w=y;for(int j=0;j<3;++j)y[j]=w[order[j]];}
        if(shape_mode>=2){
            std::array<Real,3>ell;static const int pairs[3][2]={{0,1},{0,2},{1,2}};
            for(int i=0;i<3;++i){
                if(original){for(int j=0;j<3;++j)ell[i]=ell[i]+child[i][pairs[j][0]]*child[i][pairs[j][1]]*Real(g.delta[i][j][0]*g.delta[i][j][0]);}
                else{std::array<Real,3>unit={Real(1),Real(1),Real(1)};ell[i]=common_radial_coefficients(g,g.U,alpha,child,unit,1<<i,false)[1];}
                pos(ell[i],"positive child loop variance");if(shape_mode!=5)y[i]=y[i]/ell[i];
            }
            if(shape_mode!=5){Real sum=y[0]+y[1]+y[2];for(int i=0;i<3;++i)y[i]=y[i]/sum;}
            if(mixture_branch>=0){Real E=ell[0]*y[0]+ell[1]*y[1]+ell[2]*y[2];Real ratio=root(ell[0]*ell[1]*ell[2])/(E*root(E));
                density_factor=density_factor*(mixture_branch?ratio:Real(1))/(1+ratio);}

        }
    }
    Real host[3][3];for(int e=0;e<3;++e)for(int i=0;i<3;++i)for(int j=0;j<3;++j)host[i][j]=host[i][j]+alpha[e]*Real(g.routes[e][i]*g.routes[e][j]);
    std::array<Real,2>hb={host[0][1],host[0][2]},center;
    for(int j=0;j<2;++j)for(int e=0;e<3;++e)center[j]=center[j]+Real(g.routes[e][0]*g.routes[e][j+1])/3;
    Real hostq=host[1][1]+host[2][2]-(host[1][2]+host[2][1])/2-gram_norm(hb),total;
    static const int ij[3][2]={{0,1},{0,2},{1,2}};
    for(int active=1;active<8;++active){Real q,h,de,ge,mean_de;
        if(original){Real variance[3][3];
            for(int c=0;c<3;++c)if(active&(1<<c))for(int edge=0;edge<3;++edge)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                variance[i][j]=variance[i][j]+y[c]*child[c][ij[edge][0]]*child[c][ij[edge][1]]*Real(g.delta[c][edge][i]*g.delta[c][edge][j]);
            h=variance[0][0];std::array<Real,2>cross={variance[0][1],variance[0][2]};Real energy=variance[1][1]+variance[2][2]-(variance[1][2]+variance[2][1])/2;
            q=hostq;de=energy-2*gram_dot(hb,cross)+h*gram_norm(hb);ge=h*energy-gram_norm(cross);
            mean_de=energy-2*gram_dot(center,cross)+h*gram_norm(center)+h/12;
        }else{auto u=common_radial_coefficients(g,g.U,alpha,child,y,active,false),f=common_radial_coefficients(g,g.F,alpha,child,y,active,false);
            h=u[1];q=f[0];de=f[1]-h*q;ge=f[2];if(control){auto mean=common_radial_coefficients(g,g.F,alpha,child,y,active,true);mean_de=mean[1]-h/4;}}
        total=total+Real((3-pop(active))%2?-1:1)*radial_primitive(q,h,de,ge,mean_de,control);
    }
    total=total*density_factor;if(!arb_is_finite(total.v))throw std::runtime_error("finite complete analytic radial density required");return total;
}

// Complete common-radius primitives for actual two-loop hosts.
using RadialPoly=std::array<Real,4>;
RadialPoly rp_add(const RadialPoly&a,const RadialPoly&b){RadialPoly out;for(int i=0;i<4;++i)out[i]=a[i]+b[i];return out;}
RadialPoly rp_scale(const RadialPoly&a,const Real&b){RadialPoly out;for(int i=0;i<4;++i)out[i]=a[i]*b;return out;}
RadialPoly rp_mul(const RadialPoly&a,const RadialPoly&b){RadialPoly out;for(int i=0;i<4;++i)for(int j=0;j<4-i;++j)out[i+j]=out[i+j]+a[i]*b[j];return out;}
using RadialMatrix=std::vector<std::vector<RadialPoly>>;
RadialPoly rp_minor(const RadialMatrix&m,const std::vector<int>&rows,const std::vector<int>&cols){
    int n=int(rows.size());std::vector<int>permutation(n);for(int i=0;i<n;++i)permutation[i]=i;RadialPoly sum;
    do{int inversions=0;RadialPoly term;term[0]=Real(1);
        for(int i=0;i<n;++i){for(int j=i+1;j<n;++j)if(permutation[i]>permutation[j])++inversions;term=rp_mul(term,m[rows[i]][cols[permutation[i]]]);}
        sum=rp_add(sum,rp_scale(term,Real(inversions%2?-1:1)));
    }while(std::next_permutation(permutation.begin(),permutation.end()));return sum;
}
Real radial_phi_order(const Real&t,int n){
    if(n<1||n>6)throw std::invalid_argument("bounded radial master order");pos(t,"positive master argument");
    Real pi;arb_const_pi(pi.v,precision);Real constant=sq(pi)/6;
    for(int j=1;j<n;++j)constant=constant-Real(1)/Real(j*j);
    if(arb_is_one(t.v))return constant;
    Real delta=t-1,upper=arb_cw::abs_upper(delta),half=Real(1)/2;
    if(arb_le(upper.v,half.v)){
        unsigned count=radial_series_terms(upper);Real coefficient=Real(1)/n,power=delta,sum=constant;
        for(unsigned j=1;j<=count;++j){sum=sum+Real(j%2?1:-1)*coefficient*power;
            coefficient=coefficient*Real(long(j))*Real(long(j))/(Real(long(j+1))*Real(long(n+j)));power=power*delta;}
        Real tail=coefficient*Pow::power(upper,count+1)/(1-upper);arb_cw::add_error(sum,arb_cw::abs_upper(tail));return sum;
    }
    Real dilog,arg=1-t;arb_polylog_si(dilog.v,2,arg.v,precision);Real sum=sq(pi)/6-dilog,moment=logarithm(t)/(t-1);
    for(int j=1;j<n;++j){sum=sum-t*moment/j;moment=(Real(1)/j-t*moment)/(1-t);}return sum;
}
Real proved_nonnegative(Real value,const char*label){if(arb_is_negative(value.v))throw std::runtime_error(label);arb_nonnegative_part(value.v,value.v);return value;}
std::vector<Real> radial_positive_roots(const RadialPoly&P,int degree){
    for(int j=0;j<=degree;++j)pos(P[j],"positive mass-plus-Schur coefficients");
    if(degree==2){Real disc=proved_nonnegative(sq(P[1])-4*P[0]*P[2],"quadratic real-root identity"),v=P[1]+root(disc);
        return {v/(2*P[0]),2*P[2]/v};}
    if(degree!=3)throw std::invalid_argument("proved polynomial degree required");
    Real mean=P[1]/(3*P[0]),p=proved_nonnegative((sq(P[1])-3*P[0]*P[2])/(9*sq(P[0])),"cubic real-root spread"),spread=root(p);
    if(arb_contains_zero(p.v)){
        Real low=mean-2*spread,high=mean+2*spread,hull;arb_union(hull.v,low.v,high.v,precision);pos(hull,"positive coalescent root hull");return {hull,hull,hull};
    }
    Real q=(2*P[1]*sq(P[1])-9*P[0]*P[1]*P[2]+27*sq(P[0])*P[3])/(54*P[0]*sq(P[0]));
    Real argument=q/(p*spread),minus(-1),plus(1),domain,clipped;arb_union(domain.v,minus.v,plus.v,precision);
    // The graph Gram pencil proves all three roots real; intersecting with [-1,1] retains them.
    if(!arb_intersection(clipped.v,argument.v,domain.v,precision))throw std::runtime_error("inconsistent real cubic roots");
    Real angle,pi;arb_acos(angle.v,clipped.v,precision);arb_const_pi(pi.v,precision);std::vector<Real>result;
    for(int j=0;j<3;++j){Real cosine,a=(angle+2*j*pi)/3;arb_cos(cosine.v,a.v,precision);Real value=mean+2*spread*cosine;pos(value,"positive cubic pole");result.push_back(value);}return result;
}
struct RadialMaster{
    std::vector<Real>poles,coefficients;std::vector<int>orders;Real C;bool coalescent=false;Real hull;
    explicit RadialMaster(const RadialPoly&D,int rank){
        if(rank==1){pos(D[1],"positive rank-one pole");poles={D[1]};C=logarithm(D[1])+Real(3)/2;return;}
        Real discriminant=proved_nonnegative(sq(D[1])-4*D[2],"positive generalized loop eigenvalues");Real spread=root(discriminant),hi=(D[1]+spread)/2,lo=2*D[2]/(D[1]+spread);
        pos(lo,"positive lower loop pole");pos(hi,"positive upper loop pole");poles={lo,hi};Real gap=hi-lo;
        if(arb_contains_zero(gap.v)){
            coalescent=true;arb_union(hull.v,lo.v,hi.v,precision);C=logarithm(hull);for(int j=1;j<6;++j)C=C+Real(1)/j;return;
        }
        for(int i=0;i<2;++i){Real h=poles[i],other=poles[1-i],delta=h-other;
            std::array<Real,3>co={6*h*sq(h)*sq(other)/Pow::power(delta,5),-3*h*sq(h)*other/Pow::power(delta,4),h*sq(h)/Pow::power(delta,3)};
            for(int n=1;n<=3;++n){coefficients.push_back(co[n-1]);orders.push_back(n);Real term=logarithm(h);for(int j=1;j<n;++j)term=term+Real(1)/j;C=C+co[n-1]*term;}}
    }
    Real phi(const Real&c)const{
        if(poles.size()==1)return radial_phi_order(c/poles[0],3);
        // Beta(3,3) mixture, bounded by the full positive pole hull; no point approximation.
        if(coalescent)return radial_phi_order(c/hull,6);
        Real sum;for(int i=0;i<6;++i)sum=sum+coefficients[i]*radial_phi_order(c/poles[i/3],orders[i]);return sum;
    }
};
int exact_active_loop_rank(const Graph&g,int mask){
    std::vector<std::array<int,2>>rows;for(int i=0;i<2;++i)if(mask&(1<<i))for(const auto&r:g.delta[i])if(r[0]||r[1])rows.push_back({r[0],r[1]});
    if(rows.empty())throw std::runtime_error("positive child variance required");
    for(const auto&a:rows)for(const auto&b:rows)if(a[0]*b[1]!=a[1]*b[0])return 2;return 1;
}
Real two_pole_primitive(RadialPoly D,RadialPoly F,int rank){
    if(rank==1){if(!arb_contains_zero(D[2].v)||!arb_contains_zero(F[3].v))throw std::runtime_error("exact rank-one polynomial degree");D[2]=Real(0);F[3]=Real(0);}
    D[0]=Real(1);RadialPoly P;for(int j=0;j<=rank;++j){P[j]=P[j]+D[j];P[j+1]=P[j+1]+D[j];}
    for(int j=0;j<4;++j)P[j]=P[j]+F[j];auto roots=radial_positive_roots(P,rank+1);RadialMaster master(D,rank);
    Real la=stable_log1p(F[0]),sum=sq(la)/2+master.C*la+master.phi(Real(1));
    for(const auto&h:master.poles)sum=sum+master.phi(h);for(const auto&c:roots)sum=sum-master.phi(c);return sum;
}
RadialPoly two_pole_positive(const Graph&g,const std::vector<std::pair<int,int>>&terms,const std::vector<Real>&alpha,const std::vector<std::array<Real,3>>&child,const std::array<Real,2>&shape,const Real&M,int active){
    RadialPoly result;static const int pairs[3][2]={{0,1},{0,2},{1,2}};
    for(const auto&entry:terms){Real term(entry.second);int degree=0;bool use=true;
        for(int e=0;e<6;++e)if(entry.first&(1<<e))term=term*alpha[e];
        for(int i=0;i<2;++i)for(int j=0;j<3;++j)if(entry.first&(1<<(6+3*i+j))){if(!(active&(1<<i))){use=false;break;}
            ++degree;term=term*M*shape[i]*child[i][pairs[j][0]]*child[i][pairs[j][1]];}
        if(use){if(degree>3)throw std::runtime_error("bounded radial polynomial");result[degree]=result[degree]+term;}}
    return result;
}
Real two_pole_density(const Graph&g,bool original,const std::array<Rat,11>&input){
    if(g.L!=2||g.k!=2)throw std::invalid_argument("proved two-loop host required");
    std::array<Real,11>p;for(int j=0;j<11;++j)p[j]=real(input[j]);std::array<int,3>order={0,1,2};Real total;
    do{std::array<Real,3>a;a[order[0]]=Real(1);a[order[1]]=p[0];a[order[2]]=p[0]*root(p[1]);
        std::vector<Real>alpha(6);for(int j=0;j<3;++j){alpha[g.groups[j][0]]=a[j]*p[2+j];alpha[g.groups[j][1]]=a[j]*(1-p[2+j]);}
        Real M;for(const auto&v:alpha)M=M+v;Real tropical=a[order[0]]*a[order[1]];
        std::vector<std::array<Real,3>>child={simplex(p[5],p[6]),simplex(p[7],p[8])};std::array<Real,2>shape={p[9],1-p[9]};
        Real host[4][4];for(int e=0;e<6;++e)for(int i=0;i<4;++i)for(int j=0;j<4;++j)host[i][j]=host[i][j]+alpha[e]*Real(g.routes[e][i]*g.routes[e][j]);
        Real U0=host[0][0]*host[1][1]-host[0][1]*host[1][0];pos(U0,"positive host determinant");Real sum;
        for(int active=1;active<4;++active){RadialPoly U,F;
            if(original){RadialMatrix matrix(4,std::vector<RadialPoly>(4));static const int pairs[3][2]={{0,1},{0,2},{1,2}};
                for(int i=0;i<4;++i)for(int j=0;j<4;++j){matrix[i][j][0]=host[i][j];
                    for(int c=0;c<2;++c)if(active&(1<<c))for(int e=0;e<3;++e)matrix[i][j][1]=matrix[i][j][1]+M*shape[c]*child[c][pairs[e][0]]*child[c][pairs[e][1]]*Real(g.delta[c][e][i]*g.delta[c][e][j]);}
                U=rp_minor(matrix,{0,1},{0,1});
                for(int i=0;i<2;++i)for(int j=0;j<2;++j)F=rp_add(F,rp_scale(rp_minor(matrix,{0,1,2+i},{0,1,2+j}),i==j?Real(1):-Real(1)/2));
            }else{U=two_pole_positive(g,g.U,alpha,child,shape,M,active);F=two_pole_positive(g,g.F,alpha,child,shape,M,active);}
            U=rp_scale(U,1/U0);F=rp_scale(F,1/(M*U0));int rank=exact_active_loop_rank(g,active);
            sum=sum+Real(active==3?1:-1)*two_pole_primitive(U,F,rank);}
        total=total+Pow::power(tropical/U0,3)*sum/(8*shape[0]*shape[1]);
    }while(std::next_permutation(order.begin(),order.end()));
    if(!arb_is_finite(total.v))throw std::runtime_error("finite complete two-pole radial density required");return total;
}

Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("bounded graph recipe required");auto g=load_graph(id);if(sha!=g.sha)throw std::invalid_argument("frozen graph identity required");return g;}
}
int main(int argc,char**argv){try{using namespace disjoint_triangle;
    if(argc!=2)throw std::invalid_argument("one recipe required");auto g=recipe(argv[1]);std::string mode;unsigned bits,control,reflection,transport,count;
    if(!(std::cin>>mode>>bits>>control>>reflection>>transport>>count)||(mode!="closed"&&mode!="original"&&mode!="radial-closed"&&mode!="radial-original"&&mode!="radial-dirichlet-closed"&&mode!="radial-dirichlet-original"&&mode!="radial-matched-closed"&&mode!="radial-matched-original"&&mode!="two-pole-original"&&mode!="two-pole-closed"&&mode!="radial-child-closed"&&mode!="radial-child-original"&&mode!="radial-child6-closed"&&mode!="radial-child6-original"&&mode!="radial-mixture-original"&&mode!="radial-mixture-closed")||bits<512||bits>2048||control>1||reflection<1||reflection>2||transport>2||count<1||count>256)throw std::invalid_argument("bounded native request required");precision=bits;
    std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"control\":"<<control<<",\"reflection\":"<<reflection<<",\"transport\":"<<transport<<",\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){std::array<Rat,11>p;for(auto&q:p){q=read_rational();if(!(Rat(0)<q)||!(q<Rat(1)))throw std::invalid_argument("strict open eleven cube required");}if(i)std::cout<<',';
        if(mode.rfind("radial-",0)==0){if(reflection!=1||transport!=0)throw std::invalid_argument("integrated radial mode has no radial sampling choices");
        bool original=mode.size()>=8&&mode.substr(mode.size()-8)=="original";
        int shape=mode.rfind("radial-dirichlet-",0)==0?1:(mode.rfind("radial-matched-",0)==0?2:0);
        if(mode.rfind("radial-mixture-",0)==0)interval(radial_density(g,original,bool(control),p,5,0,0)+radial_density(g,original,bool(control),p,3,0,1));
        else if(mode.rfind("radial-child6-",0)==0){Real sum;for(int j=0;j<6;++j)sum=sum+radial_density(g,original,bool(control),p,3,j);interval(sum/6);}
        else if(mode.rfind("radial-child-",0)==0)interval(radial_density(g,original,bool(control),p,3));
        else interval(radial_density(g,original,bool(control),p,shape));}
        else if(mode=="two-pole-original"||mode=="two-pole-closed"){
        if(control!=0||reflection!=1||transport!=0)throw std::invalid_argument("two-pole analytic mode: no sampled radius or conditional control");
        interval(two_pole_density(g,mode=="two-pole-original",p));}
        else interval(density(g,mode=="original",bool(control),reflection,int(transport),p));}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
