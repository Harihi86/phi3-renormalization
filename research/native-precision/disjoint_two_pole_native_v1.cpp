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
