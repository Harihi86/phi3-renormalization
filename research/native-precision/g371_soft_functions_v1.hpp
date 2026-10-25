// Keep the shifted cross propagator exact; Taylor-expand only the other factors.
#include "g371_complete_functions_v1.hpp"
std::array<Real,2> soft_host_chart(const DoubleHostGraph&g,bool positive,const V4&q,const V4&r,int excluded){
    V4 sv=sum4(q,r);Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),w=scalar4(sv,sv),X=1+w;
    pos(t,"positive soft-control q radius");pos(v,"positive soft-control r radius");pos(w,"strict intrinsic cross radius");
    std::array<Real,3>Ds={1+t,1+v,X};Real weight=power3(Ds[excluded])/(power3(Ds[0])+power3(Ds[1])+power3(Ds[2]));
    Real sqrt3=rt(Real(3));V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    auto baseline=[&](const std::array<int,4>&rr){return Real(rr[0]*rr[0])*t+Real(rr[1]*rr[1])*v+2*Real(rr[0]*rr[1])*h;};
    std::array<Real,6>den;Real d4=1;V4 hv;Real square_point,constant;std::array<V4,4>host_linear;int index=0;
    for(int j=0;j<6;++j)if(j!=2&&j!=5){den[j]=1+baseline(g.host[j]);d4=d4*den[j];V4 one=scale4(linear_vector4(g.host[j]),Real(2)/den[j]);host_linear[index++]=one;hv=sum4(hv,one);square_point=square_point+sq(scalar4(one,projection));constant=constant+external_norm4(g.host[j])/den[j];}
    V4 rv=scale4(hv,Real(-1));Real R1=scalar4(rv,projection),R2=(sq(R1)+square_point)/2-constant;
    std::array<Real,2>T0,linear,quad,con,I2;std::array<V4,2>lv,av,bv;
    for(int j=0;j<2;++j){const auto&rows=g.child[j];Real tc=baseline(rows[0]);
        if(rows[0][0]!=rows[1][0]||rows[0][1]!=rows[1][1]||rows[2][0]||rows[2][1]||!arb_equal(external_norm4(rows[2]).v,Real(1).v))throw std::invalid_argument("proved triangle pair required");
        Real a;arb_asinh(a.v,(rt(tc)/2).v,precision);Real J0=4*a/rt(tc*(tc+4)),J2=2*(1+J0)/(tc+4),ell=-(1-J0)/(2*tc),Cp=ell/2,Cf=-(J0-3*(1-J0)/tc)/12;
        T0[j]=-edge_log_full(tc)/2;I2[j]=(1-2*J0+J2)/sq(tc);av[j]=linear_vector4(rows[0]);bv[j]=linear_vector4(rows[1]);lv[j]=scale4(sum4(av[j],bv[j]),ell);
        Real qa=scalar4(av[j],projection),qb=scalar4(bv[j],projection);con[j]=Cp*(external_norm4(rows[0])+external_norm4(rows[1]))+Cf*external_norm4(rows[2]);
        linear[j]=scalar4(lv[j],projection);quad[j]=con[j]+I2[j]*(sq(qa)+qa*qb+sq(qb))/3;}
    Real P0=T0[0]*T0[1],F1=R1*P0+linear[0]*T0[1]+T0[0]*linear[1];
    Real F2=R2*P0+R1*(linear[0]*T0[1]+T0[0]*linear[1])+quad[0]*T0[1]+T0[0]*quad[1]+linear[0]*linear[1];
    Real residual,lc=2*(projection[0]+projection[1]+projection[2]+projection[3]);
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){V4 vec;for(int i=0;i<4;++i)vec[i]=Real(sign)*(Real(rr[0])*q[i]+Real(rr[1])*r[i]);vec[0]=vec[0]+Real(rr[2])-Real(rr[3])/2;vec[1]=vec[1]+sqrt3*Real(rr[3])/2;return scalar4(vec,vec);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        Real lr;for(int j=0;j<6;++j)if(j!=2&&j!=5){Real d=1+norm(g.host[j]);pos(d,"positive unexpanded host denominator");lr=lr-lp((d-den[j])/den[j]);}
        std::array<Real,2>delta;for(int j=0;j<2;++j)delta[j]=-triangle_log_full(norm(g.child[j][0]),norm(g.child[j][1]))-T0[j];
        Real difference=em1(lr)*P0+real_exp(lr)*(delta[0]*T0[1]+T0[0]*delta[1]+delta[0]*delta[1]);
        Real Y=1+norm(g.host[5]);pos(Y,"retained shifted soft propagator");residual=residual+(difference-Real(sign)*F1-F2)/(2*X*Y);}
    Real aa=w+2,root=rt(sq(w)+4),ss=aa+root,M=4*(aa+2*root)/(3*sq(ss)),mu2=aa*(aa+3*root)/(3*power3(ss));
    std::array<Real,2>along={scalar4(q,sv),scalar4(r,sv)};Real determinant=t*v-sq(h);
    if(positive){determinant=0;for(int i=0;i<4;++i)for(int j=i+1;j<4;++j)determinant=determinant+sq(q[i]*r[j]-q[j]*r[i]);}
    pos(determinant,"strict intrinsic two-vector determinant");
    V4 first={-along[0]*(aa+3*root)/(3*power3(ss)),-along[0]*(aa+3*root)/(3*power3(ss)),-along[1]*(aa+3*root)/(3*power3(ss)),-along[1]*(aa+3*root)/(3*power3(ss))};
    std::array<std::array<Real,4>,4>cov;
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int a=0;a<2;++a)for(int b=0;b<2;++b){Real gram=a==b?Real(1):Real(-1)/2,projected=mu2/4+(gram-Real(1)/4)*(M-mu2)/5;
        cov[2*i+a][2*j+b]=along[i]*along[j]*projected/w+Real(i==j?1:-1)*determinant*(gram*M-projected)/(5*w);}
    auto pair=[&](const V4&a,const V4&b){Real out;for(int i=0;i<4;++i)for(int j=0;j<4;++j)out=out+a[i]*b[j]*cov[i][j];return out;};
    Real square_mean;for(const auto&one:host_linear)square_mean=square_mean+pair(one,one);
    Real meanR2=(pair(hv,hv)+square_mean)/2-constant*M;std::array<Real,2>meanquad;
    for(int j=0;j<2;++j)meanquad[j]=con[j]*M+I2[j]*(pair(av[j],av[j])+pair(av[j],bv[j])+pair(bv[j],bv[j]))/3;
    V4 F1vector=sum4(sum4(scale4(rv,P0),scale4(lv[0],T0[1])),scale4(lv[1],T0[0]));
    Real meanF1=scalar4(F1vector,first),meanF2=meanR2*P0+pair(rv,lv[0])*T0[1]+T0[0]*pair(rv,lv[1])+meanquad[0]*T0[1]+T0[0]*meanquad[1]+pair(lv[0],lv[1]);
    Real mean=P0*(X*M-1)/sq(X)+(meanF1+meanF2)/X;
    return {weight*residual/d4,weight*mean/d4};
}
Real soft_host_value(const DoubleHostGraph&g,bool positive,bool common_scale,const std::array<Real,7>&p){
    Real t,v,radial;
    if(common_scale){Real total=sq(p[0]/(1-p[0]));t=total*p[1];v=total*(1-p[1]);radial=sq(power3(total))*sq(p[1])*sq(1-p[1])/(2*p[0]*(1-p[0]));}
    else{t=sq(p[0]/(1-p[0]));v=sq(p[1]/(1-p[1]));radial=power3(t)*power3(v)/(p[0]*(1-p[0])*p[1]*(1-p[1]));}
    Real ct,st,cp,sp,cpsi,spsi;arb_sin_cos_pi(st.v,ct.v,(2*p[4]).v,precision);arb_sin_cos_pi(sp.v,cp.v,(2*p[5]).v,precision);arb_sin_cos_pi(spsi.v,cpsi.v,p[6].v,precision);
    Real z=p[2],Z=p[3];V4 e={rt(z)*ct,rt(z)*st,rt(1-z),Real(0)},n1={-st,ct,Real(0),Real(0)},n2={rt(1-z)*ct,rt(1-z)*st,-rt(z),Real(0)},n3={Real(0),Real(0),Real(0),Real(1)};
    V4 n=sum4(scale4(sum4(scale4(n1,cp),scale4(n2,sp)),rt(Z)),scale4(n3,rt(1-Z)));
    V4 q=scale4(e,rt(t)),r=scale4(sum4(scale4(e,cpsi),scale4(n,spsi)),rt(v));
    auto a=soft_host_chart(g,positive,q,r,2),b=soft_host_chart(g,positive,q,sum4(r,q,-1),1),c=soft_host_chart(g,positive,sum4(r,q,-1),q,0);
    Real residual=a[0]+b[0]+c[0],mean=a[1]+b[1]+c[1];Real orientation=3*(1-z)*rt(1-Z),intrinsic=Real(8)/3*sq(sq(spsi));
    Real answer=radial*intrinsic*(orientation*residual+mean);if(!arb_is_finite(answer.v))throw std::runtime_error("finite soft-propagator control required");return answer;
}
