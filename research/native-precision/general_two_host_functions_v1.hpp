// Inside nested_closed after the frozen two-host coordinate helpers.
#include "general_two_host_data.hpp"
Real inverse_triangle_general(const Real&a,const Real&b,const Real&c){
    pos(a,"positive general triangle a");pos(b,"positive general triangle b");pos(c,"positive general triangle c");
    Real alpha=sq(a-b)-2*c*(a+b)+sq(c),beta=2*c*(a+b-c),root=rt(c*(c+4));pos(-alpha,"strict general triangle Gram");
    Real kappa=rt(-alpha),z0=beta/(2*root),q1=rt(sq(a-b)+4*c),delta=(4*sq(root)*alpha-sq(beta))/(2*root*(2*root*q1+2*sq(root)+beta));
    pos(-delta,"decreasing general conic chart");Z pole;arb_set(acb_imagref(pole.v),kappa.v);
    Z kernel=Z(delta)/(Z(z0)-pole),lk;acb_log1p(lk.v,kernel.v,precision);
    Z result=Z(2*real_log((root+c)/(root-c)))*lk;
    std::array<Real,2>linear={a-b-c,-a+b-c},square={b,a};
    for(int j=0;j<2;++j)for(long sign:{1L,-1L}){
        Real leading=Real(sign)*c-root,middle=beta-2*root*Real(sign)*linear[j],disc=4*c*rt(square[j]*(square[j]+4));
        for(long side:{1L,-1L}){Real zero=(-middle+Real(side)*disc)/(2*leading),factor=delta/(z0-zero);result=result+Z(sign)*pair_log_integral(factor,kernel);}}
    Real im;arb_set(im.v,acb_imagref(result.v));return -2*im/kappa;
}
Real triangle_log_general(const Real&a,const Real&b,const Real&c){
    if(arb_is_zero(c.v)){if(!arb_equal(a.v,b.v))throw std::invalid_argument("exact equal degenerate triangle invariants required");return edge_log_full(a)/2;}
    if(arb_equal(c.v,Real(1).v))return triangle_log_full(a,b);
    Real cross=(a+b-c)/2,det=a*b-sq(cross);pos(det,"strict general Euclidean Gram");
    Real cx=b*(a-cross)/(2*det),cy=a*(b-cross)/(2*det),height=1+(a*cx+b*cy)/2;
    return (cx*edge_log_full(b)+cy*edge_log_full(a)+(1-cx-cy)*edge_log_full(c))/2-Real(1)/2+height*inverse_triangle_general(a,b,c);
}
Real general_host_chart(const DoubleHostGraph&g,bool positive,const V4&q,const V4&r,int excluded){
    Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),sum=t+v+2*h;
    pos(t,"positive first host radius");pos(v,"positive second host radius");pos(sum,"positive third host radius");
    std::array<Real,3>Ds={1+t,1+v,1+sum};Real weight=power3(Ds[excluded])/(power3(Ds[0])+power3(Ds[1])+power3(Ds[2]));
    Real sqrt3=rt(Real(3));V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    auto baseline=[&](const std::array<int,4>&rr){return Real(rr[0]*rr[0])*t+Real(rr[1]*rr[1])*v+2*Real(rr[0]*rr[1])*h;};
    std::array<Real,6>den;Real denominator=1;
    for(int j=0;j<6;++j){den[j]=1+baseline(g.host[j]);pos(den[j],"positive baseline host denominator");denominator=denominator*den[j];}
    std::array<Real,2>T0;
    for(int j=0;j<2;++j){const auto&rows=g.child[j];T0[j]=-triangle_log_general(baseline(rows[0]),baseline(rows[1]),baseline(rows[2]));}
    Real total,P0=T0[0]*T0[1];
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){V4 vec;for(int i=0;i<4;++i)vec[i]=Real(sign)*(Real(rr[0])*q[i]+Real(rr[1])*r[i]);
                vec[0]=vec[0]+Real(rr[2])-Real(rr[3])/2;vec[1]=vec[1]+sqrt3*Real(rr[3])/2;return scalar4(vec,vec);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        Real lr;for(int j=0;j<6;++j){Real d=1+norm(g.host[j]);pos(d,"positive physical host denominator");lr=lr-lp((d-den[j])/den[j]);}
        std::array<Real,2>delta;for(int j=0;j<2;++j){const auto&rows=g.child[j];delta[j]=-triangle_log_general(norm(rows[0]),norm(rows[1]),norm(rows[2]))-T0[j];}
        total=total+(em1(lr)*P0+real_exp(lr)*(delta[0]*T0[1]+T0[0]*delta[1]+delta[0]*delta[1]))/2;
    }
    return weight*total/denominator;
}
Real general_host_value(const DoubleHostGraph&g,bool positive,bool common_scale,const std::array<Real,7>&p){
    Real t,v,measure;
    if(common_scale){Real total=sq(p[0]/(1-p[0]));t=total*p[1];v=total*(1-p[1]);measure=4*sq(power3(total))*sq(p[1])*sq(1-p[1])/(p[0]*(1-p[0]));}
    else{t=sq(p[0]/(1-p[0]));v=sq(p[1]/(1-p[1]));measure=8*power3(t)*power3(v)/(p[0]*(1-p[0])*p[1]*(1-p[1]));}
    Real c1,s1,c2,s2,cw,sw;arb_sin_cos_pi(s1.v,c1.v,(2*p[4]).v,precision);arb_sin_cos_pi(s2.v,c2.v,(2*p[5]).v,precision);arb_sin_cos_pi(sw.v,cw.v,p[6].v,precision);
    V4 u={rt(t*p[2])*c1,rt(t*p[2])*s1,rt(t*(1-p[2])),Real(0)};
    V4 vv={rt(v*p[3])*c2,rt(v*p[3])*s2,rt(v*(1-p[3]))*cw,rt(v*(1-p[3]))*sw};
    Real value=general_host_chart(g,positive,u,vv,2)+general_host_chart(g,positive,u,sum4(vv,u,-1),1)+general_host_chart(g,positive,sum4(vv,u,-1),u,0);
    value=value*measure*(1-p[2])*(1-p[3])*sq(sw);if(!arb_is_finite(value.v))throw std::runtime_error("finite general two-host forest required");return value;
}
DoubleHostGraph general_host_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original graph recipe");auto g=load_general_host(id);if(sha!=g.sha)throw std::invalid_argument("exact original two-host identity required");return g;}
