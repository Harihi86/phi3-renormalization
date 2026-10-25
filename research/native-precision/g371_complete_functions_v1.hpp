// Two complete triangle children. Included in nested_closed after scalar helpers.
struct DoubleHostGraph{std::string id,sha;std::array<std::array<int,4>,6>host;std::array<std::array<std::array<int,4>,3>,2>child;};
#include "g371_complete_data.hpp"
using V4=std::array<Real,4>;
Real scalar4(const V4&a,const V4&b){Real value;for(int i=0;i<4;++i)value=value+a[i]*b[i];return value;}
V4 sum4(const V4&a,const V4&b,long sign=1){V4 v;for(int i=0;i<4;++i)v[i]=a[i]+Real(sign)*b[i];return v;}
V4 scale4(const V4&a,const Real&b){V4 v;for(int i=0;i<4;++i)v[i]=a[i]*b;return v;}
Real external_norm4(const std::array<int,4>&r){return Real(r[2]*r[2]+r[3]*r[3]-r[2]*r[3]);}
V4 linear_vector4(const std::array<int,4>&r){return {Real(r[0]*r[2]),Real(r[0]*r[3]),Real(r[1]*r[2]),Real(r[1]*r[3])};}
Real rotational_pair(const V4&a,const V4&b,const Real&t,const Real&v,const Real&h){
    auto g=[&](int i,int j){return a[i]*b[j]+a[i+1]*b[j+1]-(a[i]*b[j+1]+a[i+1]*b[j])/2;};
    return (t*g(0,0)+v*g(2,2)+h*(g(0,2)+g(2,0)))/6;
}
Real double_host_chart(const DoubleHostGraph&g,bool positive,bool control,const V4&q,const V4&r,int excluded){
    Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),sum=t+v+2*h;
    pos(t,"positive first host radius");pos(v,"positive second host radius");pos(sum,"positive third host radius");
    std::array<Real,3>Ds={1+t,1+v,1+sum};Real weight=power3(Ds[excluded])/(power3(Ds[0])+power3(Ds[1])+power3(Ds[2]));
    Real sqrt3=rt(Real(3));V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    auto baseline=[&](const std::array<int,4>&rr){return Real(rr[0]*rr[0])*t+Real(rr[1]*rr[1])*v+2*Real(rr[0]*rr[1])*h;};
    std::array<Real,6>den;Real denominator=1;
    for(int j=0;j<6;++j){den[j]=1+baseline(g.host[j]);pos(den[j],"positive baseline host denominator");denominator=denominator*den[j];}
    std::array<Real,2>T0,ell,Cp,Cf,I2,tc;
    for(int j=0;j<2;++j){const auto&rows=g.child[j];tc[j]=baseline(rows[0]);
        if(rows[0][0]!=rows[1][0]||rows[0][1]!=rows[1][1]||rows[2][0]||rows[2][1]||!arb_equal(external_norm4(rows[2]).v,Real(1).v))throw std::invalid_argument("proved triangle loop and fixed invariant required");
        Real a;arb_asinh(a.v,(rt(tc[j])/2).v,precision);Real J0=4*a/rt(tc[j]*(tc[j]+4)),J2=2*(1+J0)/(tc[j]+4);
        T0[j]=-edge_log_full(tc[j])/2;ell[j]=-(1-J0)/(2*tc[j]);Cp[j]=ell[j]/2;
        Cf[j]=-(J0-3*(1-J0)/tc[j])/12;I2[j]=(1-2*J0+J2)/sq(tc[j]);}
    Real total,P0=T0[0]*T0[1];
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){V4 vec;for(int i=0;i<4;++i)vec[i]=Real(sign)*(Real(rr[0])*q[i]+Real(rr[1])*r[i]);
                vec[0]=vec[0]+Real(rr[2])-Real(rr[3])/2;vec[1]=vec[1]+sqrt3*Real(rr[3])/2;return scalar4(vec,vec);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        Real logratio;for(int j=0;j<6;++j){Real d=1+norm(g.host[j]);pos(d,"positive physical host denominator");logratio=logratio-lp((d-den[j])/den[j]);}
        std::array<Real,2>delta;for(int j=0;j<2;++j)delta[j]=-triangle_log_full(norm(g.child[j][0]),norm(g.child[j][1]))-T0[j];
        total=total+(em1(logratio)*P0+real_exp(logratio)*(delta[0]*T0[1]+T0[0]*delta[1]+delta[0]*delta[1]))/2;
    }
    if(control){
        V4 hv;Real square_point,square_mean,constant;
        for(int j=0;j<6;++j){V4 one=scale4(linear_vector4(g.host[j]),Real(2)/den[j]);hv=sum4(hv,one);
            square_point=square_point+sq(scalar4(one,projection));square_mean=square_mean+rotational_pair(one,one,t,v,h);constant=constant+external_norm4(g.host[j])/den[j];}
        V4 rv=scale4(hv,Real(-1));Real R1=scalar4(rv,projection),R2=(sq(R1)+square_point)/2-constant;
        Real meanR2=(rotational_pair(hv,hv,t,v,h)+square_mean)/2-constant;
        std::array<V4,2>lv;std::array<Real,2>linear,quad,meanquad;
        for(int j=0;j<2;++j){const auto&rows=g.child[j];V4 av=linear_vector4(rows[0]),bv=linear_vector4(rows[1]);lv[j]=scale4(sum4(av,bv),ell[j]);
            Real qa=scalar4(av,projection),qb=scalar4(bv,projection),con=Cp[j]*(external_norm4(rows[0])+external_norm4(rows[1]))+Cf[j]*external_norm4(rows[2]);
            linear[j]=scalar4(lv[j],projection);quad[j]=con+I2[j]*(sq(qa)+qa*qb+sq(qb))/3;
            meanquad[j]=con+I2[j]*(rotational_pair(av,av,t,v,h)+rotational_pair(av,bv,t,v,h)+rotational_pair(bv,bv,t,v,h))/3;}
        Real point=R2*P0+R1*(linear[0]*T0[1]+T0[0]*linear[1])+quad[0]*T0[1]+T0[0]*quad[1]+linear[0]*linear[1];
        Real mean=meanR2*P0+rotational_pair(rv,lv[0],t,v,h)*T0[1]+T0[0]*rotational_pair(rv,lv[1],t,v,h)+meanquad[0]*T0[1]+T0[0]*meanquad[1]+rotational_pair(lv[0],lv[1],t,v,h);
        total=total+mean-point;
    }
    return weight*total/denominator;
}
Real double_host_value(const DoubleHostGraph&g,bool positive,bool common_scale,bool control,const std::array<Real,7>&p){
    Real t,v,measure;
    if(common_scale){Real total=sq(p[0]/(1-p[0]));t=total*p[1];v=total*(1-p[1]);measure=4*sq(power3(total))*sq(p[1])*sq(1-p[1])/(p[0]*(1-p[0]));}
    else{t=sq(p[0]/(1-p[0]));v=sq(p[1]/(1-p[1]));measure=8*power3(t)*power3(v)/(p[0]*(1-p[0])*p[1]*(1-p[1]));}
    Real c1,s1,c2,s2,cw,sw;arb_sin_cos_pi(s1.v,c1.v,(2*p[4]).v,precision);arb_sin_cos_pi(s2.v,c2.v,(2*p[5]).v,precision);arb_sin_cos_pi(sw.v,cw.v,p[6].v,precision);
    V4 u={rt(t*p[2])*c1,rt(t*p[2])*s1,rt(t*(1-p[2])),Real(0)};
    V4 vv={rt(v*p[3])*c2,rt(v*p[3])*s2,rt(v*(1-p[3]))*cw,rt(v*(1-p[3]))*sw};
    Real value=double_host_chart(g,positive,control,u,vv,2)+double_host_chart(g,positive,control,u,sum4(vv,u,-1),1)+double_host_chart(g,positive,control,sum4(vv,u,-1),u,0);
    value=value*measure*(1-p[2])*(1-p[3])*sq(sw);
    if(!arb_is_finite(value.v))throw std::runtime_error("finite complete seven-cube forest required");return value;
}
DoubleHostGraph double_host_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original graph recipe");auto g=make_g371_complete();if(id!=g.id||sha!=g.sha)throw std::invalid_argument("exact original G371 identity required");return g;}
