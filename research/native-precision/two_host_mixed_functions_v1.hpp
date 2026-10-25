// Complete bubble plus triangle after exact original-forest elimination.
struct TwoHostMixedGraph{std::string id,sha;std::array<std::array<int,4>,6>host;std::array<std::array<int,4>,3>triangle;std::array<int,4>bubble;};
#include "two_host_mixed_data.hpp"
std::array<Real,3> two_host_bubble_jets(const Real&t){
    pos(t,"positive bubble baseline");Real E=edge_log_full(t),a;arb_asinh(a.v,(rt(t)/2).v,precision);
    Real j=4*a/rt(t*(t+4)),B=(t+4)*E/6-t/9,B1=E/6+(t+4)*(1-j)/(6*t)-Real(1)/9,B2=Real(1)/(6*t)-(1-j)/sq(t),D=1+t;
    return {B/D,B1/D-B/sq(D),B2/D-2*B1/sq(D)+2*B/power3(D)};
}
Real two_host_bubble_value(const Real&t){pos(t,"positive shifted bubble");return ((t+4)*edge_log_full(t)/6-t/9)/(1+t);}
std::array<Real,2> two_host_mixed_chart(const TwoHostMixedGraph&g,bool positive,const V4&q,const V4&r,int excluded){
    Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),w=scalar4(sum4(q,r),sum4(q,r));
    pos(t,"positive mixed q");pos(v,"positive mixed r");pos(w,"positive mixed sum");
    std::array<Real,3>ds={1+t,1+v,1+w};Real weight=power3(ds[excluded])/(power3(ds[0])+power3(ds[1])+power3(ds[2]));
    Real sqrt3=rt(Real(3));V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    auto baseline=[&](const std::array<int,4>&rr){return Real(rr[0]*rr[0])*t+Real(rr[1]*rr[1])*v+2*Real(rr[0]*rr[1])*h;};
    std::array<Real,6>den;Real denominator=1;QuadraticJet product;product.value=1;
    for(int k=0;k<6;++k){den[k]=1+baseline(g.host[k]);pos(den[k],"positive baseline mixed host");denominator=denominator*den[k];
        QuadraticJet one;one.value=1;one.linear=scale4(linear_vector4(g.host[k]),Real(-2)/den[k]);one.constant=-external_norm4(g.host[k])/den[k];
        for(int i=0;i<4;++i)for(int j=0;j<4;++j)one.square[i][j]=one.linear[i]*one.linear[j];product=multiply_quadratic(product,one);}
    auto bj=two_host_bubble_jets(baseline(g.bubble));QuadraticJet bubble;bubble.value=bj[0];
    V4 bl=scale4(linear_vector4(g.bubble),Real(2));bubble.linear=scale4(bl,bj[1]);bubble.constant=bj[1]*external_norm4(g.bubble);
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)bubble.square[i][j]=bj[2]*bl[i]*bl[j]/2;
    product=multiply_quadratic(product,bubble);
    const auto&rows=g.triangle;QuadraticJet child;std::array<Real,3>base={baseline(rows[0]),baseline(rows[1]),baseline(rows[2])};
    if(rows[2][0]==0&&rows[2][1]==0){
        if(rows[0][0]!=rows[1][0]||rows[0][1]!=rows[1][1]||!arb_equal(external_norm4(rows[2]).v,Real(1).v))throw std::invalid_argument("exact unit mixed triangle required");
        Real tc=base[0],as;arb_asinh(as.v,(rt(tc)/2).v,precision);Real J0=4*as/rt(tc*(tc+4)),J2=2*(1+J0)/(tc+4),ell=-(1-J0)/(2*tc),Cp=ell/2,Cf=-(J0-3*(1-J0)/tc)/12,I2=(1-2*J0+J2)/sq(tc);
        V4 av=linear_vector4(rows[0]),bv=linear_vector4(rows[1]);child.value=-edge_log_full(tc)/2;child.linear=scale4(sum4(av,bv),ell);
        child.constant=Cp*(external_norm4(rows[0])+external_norm4(rows[1]))+Cf*external_norm4(rows[2]);
        for(int i=0;i<4;++i)for(int j=0;j<4;++j)child.square[i][j]=I2*(av[i]*av[j]+bv[i]*bv[j]+(av[i]*bv[j]+bv[i]*av[j])/2)/3;
    }else{
        auto jet=general_triangle_jet(base[0],base[1],base[2]);child.value=-jet.value;
        std::array<V4,3>vectors={linear_vector4(rows[0]),linear_vector4(rows[1]),linear_vector4(rows[2])};
        for(int a=0;a<3;++a){child.linear=sum4(child.linear,scale4(vectors[a],-2*jet.gradient[a]));child.constant=child.constant-jet.gradient[a]*external_norm4(rows[a]);
            for(int b=0;b<3;++b)for(int i=0;i<4;++i)for(int j=0;j<4;++j)child.square[i][j]=child.square[i][j]-2*jet.hessian[a][b]*vectors[a][i]*vectors[b][j];}
    }
    product=multiply_quadratic(product,child);Real total,P0=bubble.value*child.value;
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){V4 vec;for(int i=0;i<4;++i)vec[i]=Real(sign)*(Real(rr[0])*q[i]+Real(rr[1])*r[i]);vec[0]=vec[0]+Real(rr[2])-Real(rr[3])/2;vec[1]=vec[1]+sqrt3*Real(rr[3])/2;return scalar4(vec,vec);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        Real lr;for(int k=0;k<6;++k){Real dd=1+norm(g.host[k]);pos(dd,"positive shifted mixed host");lr=lr-lp((dd-den[k])/den[k]);}
        Real db=two_host_bubble_value(norm(g.bubble))-bubble.value;
        Real dt=-triangle_log_general(norm(rows[0]),norm(rows[1]),norm(rows[2]))-child.value;
        total=total+(em1(lr)*P0+real_exp(lr)*(db*child.value+bubble.value*dt+db*dt))/2;
    }
    Real point=product.constant,mean=product.constant;
    for(int i=0;i<4;++i)for(int j=0;j<4;++j){Real loop=(i/2==0&&j/2==0)?t:(i/2==1&&j/2==1)?v:h;Real external=i%2==j%2?Real(1):Real(-1)/2;
        point=point+product.square[i][j]*projection[i]*projection[j];mean=mean+product.square[i][j]*loop*external/6;}
    return {weight*(total-point)/denominator,weight*mean/denominator};
}
TwoHostMixedGraph two_host_mixed_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one mixed graph recipe required");auto g=load_two_host_mixed(id);if(sha!=g.sha)throw std::invalid_argument("exact original mixed identity required");return g;}
