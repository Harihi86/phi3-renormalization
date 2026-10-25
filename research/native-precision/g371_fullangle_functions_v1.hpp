// Exact positive-measure marginal of the complete, unpartitioned G371 kernel.
#include "g371_complete_functions_v1.hpp"
Real fullangle_host_value(const DoubleHostGraph&g,bool positive,bool common_scale,bool control,const std::array<Real,6>&p){
    Real t,v,measure;
    if(common_scale){Real total=sq(p[0]/(1-p[0]));t=total*p[1];v=total*(1-p[1]);measure=2*sq(power3(total))*sq(p[1])*sq(1-p[1])/(p[0]*(1-p[0]));}
    else{t=sq(p[0]/(1-p[0]));v=sq(p[1]/(1-p[1]));measure=4*power3(t)*power3(v)/(p[0]*(1-p[0])*p[1]*(1-p[1]));}
    Real c1,s1,c2,s2;arb_sin_cos_pi(s1.v,c1.v,(2*p[4]).v,precision);arb_sin_cos_pi(s2.v,c2.v,(2*p[5]).v,precision);
    Real qx=rt(t*p[2])*c1,qy=rt(t*p[2])*s1,rx=rt(v*p[3])*c2,ry=rt(v*p[3])*s2;
    Real qp=rt(t*(1-p[2])),rp=rt(v*(1-p[3])),sqrt3=rt(Real(3));
    V4 projection={qx,(-qx+sqrt3*qy)/2,rx,(-rx+sqrt3*ry)/2};
    Real B=2*qp*rp,A=1+t+v+2*(qx*rx+qy*ry),am=A-B,ap=A+B;
    if(positive){am=1+sq(qx+rx)+sq(qy+ry)+sq(qp-rp);ap=1+sq(qx+rx)+sq(qy+ry)+sq(qp+rp);A=(am+ap)/2;}
    pos(am,"strict massive baseline angle gap");pos(ap,"strict positive baseline angle sum");
    Real ra=rt(am*ap),m1=2/(A+ra),m2=2/(ra*(A+ra)),m3=1/power3(ra),m4=A/(ra*sq(sq(ra)));
    auto baseline=[&](const std::array<int,4>&rr){if(rr[0]*rr[1])throw std::invalid_argument("angle-free row required");return Real(rr[0]*rr[0])*t+Real(rr[1]*rr[1])*v;};
    std::array<Real,6>den;Real d4=1;
    for(int j=0;j<6;++j)if(j!=2&&j!=5){den[j]=1+baseline(g.host[j]);d4=d4*den[j];}
    std::array<Real,2>T0,ell,Cp,Cf,I2;
    for(int j=0;j<2;++j){const auto&rows=g.child[j];Real tc=baseline(rows[0]);
        if(rows[0][0]!=rows[1][0]||rows[0][1]!=rows[1][1]||rows[2][0]||rows[2][1]||!arb_equal(external_norm4(rows[2]).v,Real(1).v))throw std::invalid_argument("proved triangle pair required");
        T0[j]=-edge_log_full(tc)/2;
        if(control){Real a;arb_asinh(a.v,(rt(tc)/2).v,precision);Real J0=4*a/rt(tc*(tc+4)),J2=2*(1+J0)/(tc+4);
            ell[j]=-(1-J0)/(2*tc);Cp[j]=ell[j]/2;Cf[j]=-(J0-3*(1-J0)/tc)/12;I2[j]=(1-2*J0+J2)/sq(tc);}}
    Real global_mean;Real P0=T0[0]*T0[1],total,lc=2*(projection[0]+projection[1]+projection[2]+projection[3]);
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){Real x=Real(sign)*(Real(rr[0])*qx+Real(rr[1])*rx)+Real(rr[2])-Real(rr[3])/2;
                Real y=Real(sign)*(Real(rr[0])*qy+Real(rr[1])*ry)+sqrt3*Real(rr[3])/2;
                return sq(x)+sq(y)+Real(rr[0]*rr[0])*sq(qp)+Real(rr[1]*rr[1])*sq(rp);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        Real C=A+Real(sign)*lc+1,cm=C-B,cp=C+B;
        if(positive){Real x=Real(sign)*(qx+rx)+Real(1)/2,y=Real(sign)*(qy+ry)+sqrt3/2;
            cm=1+sq(x)+sq(y)+sq(qp-rp);cp=1+sq(x)+sq(y)+sq(qp+rp);C=(cm+cp)/2;}
        pos(cm,"strict massive shifted angle gap");pos(cp,"strict positive shifted angle sum");
        Real rc=rt(cm*cp),K=2*(1+(A+C)/(ra+rc))/((A+ra)*(C+rc));
        Real lr;for(int j=0;j<6;++j)if(j!=2&&j!=5){Real d=1+norm(g.host[j]);pos(d,"positive other host denominator");lr=lr-lp((d-den[j])/den[j]);}
        std::array<Real,2>delta;for(int j=0;j<2;++j)delta[j]=-triangle_log_full(norm(g.child[j][0]),norm(g.child[j][1]))-T0[j];
        Real difference=em1(lr)*P0+real_exp(lr)*(delta[0]*T0[1]+T0[0]*delta[1]+delta[0]*delta[1]);
        total=total+(P0*(K-m2)+K*difference)/2;
    }
    if(control){
        V4 hv;Real square_point,constant;std::array<V4,4>host_linear;int index=0;
        for(int j=0;j<6;++j)if(j!=2&&j!=5){V4 one=scale4(linear_vector4(g.host[j]),Real(2)/den[j]);host_linear[index++]=one;hv=sum4(hv,one);
            square_point=square_point+sq(scalar4(one,projection));constant=constant+external_norm4(g.host[j])/den[j];}
        V4 rv=scale4(hv,Real(-1));Real R1=scalar4(rv,projection),R2=(sq(R1)+square_point)/2-constant;
        std::array<V4,2>lv,av,bv;std::array<Real,2>linear,quad,con;
        for(int j=0;j<2;++j){const auto&rows=g.child[j];av[j]=linear_vector4(rows[0]);bv[j]=linear_vector4(rows[1]);lv[j]=scale4(sum4(av[j],bv[j]),ell[j]);
            Real qa=scalar4(av[j],projection),qb=scalar4(bv[j],projection);
            con[j]=Cp[j]*(external_norm4(rows[0])+external_norm4(rows[1]))+Cf[j]*external_norm4(rows[2]);
            linear[j]=scalar4(lv[j],projection);quad[j]=con[j]+I2[j]*(sq(qa)+qa*qb+sq(qb))/3;}
        Real F1=R1*P0+linear[0]*T0[1]+T0[0]*linear[1];
        Real F2=R2*P0+R1*(linear[0]*T0[1]+T0[0]*linear[1])+quad[0]*T0[1]+T0[0]*quad[1]+linear[0]*linear[1];
        V4 f1vector=sum4(sum4(scale4(rv,P0),scale4(lv[0],T0[1])),scale4(lv[1],T0[0]));
        V4 cross={Real(2),Real(2),Real(2),Real(2)};
        auto meanF2=[&](const Real&h){
            Real square_mean;for(const auto&one:host_linear)square_mean=square_mean+rotational_pair(one,one,t,v,h);
            Real meanR2=(rotational_pair(hv,hv,t,v,h)+square_mean)/2-constant;
            std::array<Real,2>meanquad;for(int j=0;j<2;++j)meanquad[j]=con[j]+I2[j]*(rotational_pair(av[j],av[j],t,v,h)+rotational_pair(av[j],bv[j],t,v,h)+rotational_pair(bv[j],bv[j],t,v,h))/3;
            return meanR2*P0+rotational_pair(rv,lv[0],t,v,h)*T0[1]+T0[0]*rotational_pair(rv,lv[1],t,v,h)+meanquad[0]*T0[1]+T0[0]*meanquad[1]+rotational_pair(lv[0],lv[1],t,v,h);};
        Real aF=meanF2(Real(0)),bF=meanF2(Real(1))-aF;
        Real aL=rotational_pair(cross,f1vector,t,v,Real(0)),bL=rotational_pair(cross,f1vector,t,v,Real(1))-aL;
        Real point=m2*F2-m3*(lc*F1+P0)+m4*sq(lc)*P0;
        Real base=1+t+v,fullb=2*rt(t*v),gap=1+sq(rt(t)-rt(v));
        Real fullr=rt(gap*(base+fullb)),S=base+fullr;
        Real L1=4*(base+2*fullr)/(3*sq(S)),L2=4/sq(S),L3=4/(fullr*sq(S)),L4=4*(base+2*fullr)/(3*power3(fullr)*sq(S));
        global_mean=aF*L2+bF*(L1-base*L2)/2-aL*L3-bL*(L2-base*L3)/2-P0*L3+Real(2)/3*P0*(L3-L4);
        total=total-point;
    }
    Real answer=measure*((1-p[2])*(1-p[3])*total+global_mean/4)/d4;
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite six-cube complete angular marginal required");return answer;
}
