// Included inside nested_closed. Reuses only proved scalar triangle helpers.
struct TripleGraph{std::string id,sha;std::array<std::array<int,3>,3>host;std::array<std::array<std::array<int,3>,3>,3>child;};
#include "g521_complete_data.hpp"
Real triple_value(const TripleGraph&g,bool positive,bool control,const std::array<Real,3>&p){
    Real s=p[0],z=p[1],theta=p[2],t=sq(s/(1-s)),A=1+t,co,si;
    arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real rad=rt(t*z),sqrt3=rt(Real(3)),T0=-edge_log_full(t)/2,total;
    for(long sign:{-1L,1L}){
        Real qx=Real(sign)*rad*co,qy=Real(sign)*rad*si,q1=qx,q2=(-qx+sqrt3*qy)/2;
        auto norm=[&](const std::array<int,3>&r){
            if(positive)return sq(Real(r[0])*qx+Real(r[1])-Real(r[2])/2)+sq(Real(r[0])*qy+sqrt3*Real(r[2])/2)+Real(r[0]*r[0])*t*(1-z);
            return Real(r[0]*r[0])*t+2*Real(r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real lr;for(const auto&r:g.host){Real prop=1+norm(r);pos(prop,"positive original host propagator");lr=lr-lp((prop-A)/A);}
        std::array<Real,3>d;
        for(int i=0;i<3;++i){const auto&child=g.child[i];
            for(int j=0;j<3;++j)if(child[j][0]*child[j][0]!=(j<2?1:0))throw std::invalid_argument("proved child loop partition required");
            const auto&fixed=child[2];if(fixed[1]*fixed[1]+fixed[2]*fixed[2]-fixed[1]*fixed[2]!=1)throw std::invalid_argument("unit original child invariant required");
            d[i]=-triangle_log_full(norm(child[0]),norm(child[1]))-T0;
        }
        Real difference=sq(T0)*(d[0]+d[1]+d[2])+T0*(d[0]*d[1]+d[0]*d[2]+d[1]*d[2])+d[0]*d[1]*d[2];
        total=total+(em1(lr)*power3(T0)+real_exp(lr)*difference)/2;
    }
    Real common=2/(s*(1-s))*power3(t/A),answer=common*(1-z)*total;
    if(control){
        std::array<Real,2>q={rad*co,rad*(-co+sqrt3*si)/2},hv;Real hsquare,hconstant,hpoint;
        for(const auto&r:g.host){std::array<Real,2>one={Real(2*r[0]*r[1]),Real(2*r[0]*r[2])};
            hsquare=hsquare+norm_gram(one);hconstant=hconstant+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);
            hpoint=hpoint+sq(q[0]*one[0]+q[1]*one[1]);for(int j=0;j<2;++j)hv[j]=hv[j]+one[j];}
        std::array<Real,2>rv={-hv[0]/A,-hv[1]/A};Real R1=q[0]*rv[0]+q[1]*rv[1];
        Real R2=(sq(q[0]*hv[0]+q[1]*hv[1])+hpoint)/(2*sq(A))-hconstant/A;
        Real meanR2=t*(norm_gram(hv)+hsquare)/(12*sq(A))-hconstant/A;
        Real h;arb_asinh(h.v,(rt(t)/2).v,precision);Real J0=4*h/rt(t*(t+4)),J2=2*(1+J0)/(t+4);
        Real U1=-(1-J0)/(2*t),Cpair=U1/2,Cfixed=-(J0-3*(1-J0)/t)/12,I2=(1-2*J0+J2)/sq(t);
        std::array<std::array<Real,2>,3>vectors;std::array<Real,3>linear,quadratic,meanQuadratic;
        for(int i=0;i<3;++i){const auto&rows=g.child[i];std::array<Real,2>a,b;
            for(int j=0;j<2;++j){a[j]=Real(rows[0][0]*rows[0][j+1]);b[j]=Real(rows[1][0]*rows[1][j+1]);vectors[i][j]=a[j]+b[j];}
            Real qa=q[0]*a[0]+q[1]*a[1],qb=q[0]*b[0]+q[1]*b[1],constant;
            for(int j=0;j<3;++j){const auto&r=rows[j];constant=constant+(j<2?Cpair:Cfixed)*Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);}
            linear[i]=U1*(qa+qb);quadratic[i]=constant+I2*(sq(qa)+qa*qb+sq(qb))/3;
            meanQuadratic[i]=constant+I2*t*(norm_gram(a)+dot_gram(a,b)+norm_gram(b))/18;
        }
        Real point=R2*power3(T0),mean=meanR2*power3(T0);
        for(int i=0;i<3;++i){point=point+sq(T0)*(R1*linear[i]+quadratic[i]);mean=mean+sq(T0)*(U1*t*dot_gram(rv,vectors[i])/6+meanQuadratic[i]);
            for(int j=i+1;j<3;++j){point=point+T0*linear[i]*linear[j];mean=mean+T0*sq(U1)*t*dot_gram(vectors[i],vectors[j])/6;}}
        answer=answer+common*(mean/2-(1-z)*point);
    }
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite complete original three-cube forest required");return answer;
}
TripleGraph triple_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original graph recipe");auto g=make_g521_complete();if(id!=g.id||sha!=g.sha)throw std::invalid_argument("exact original G521 identity required");return g;}
