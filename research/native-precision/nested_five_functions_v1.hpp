// Included inside nested_closed, after the verified six-cube helpers.
Real edge_log_full(const Real&t){return 2*rt((t+4)/t)*([](const Real&x){Real r;arb_asinh(r.v,x.v,precision);return r;})(rt(t)/2)-2;}
Z pair_log_integral(const Real&a,const Z&b){
    pos(1+a,"positive real logarithm ratio");
    if(arb_contains_zero(acb_imagref(b.v)))throw std::runtime_error("strict off-real kernel required");
    if(arb_is_zero(a.v))return Z(0);
    Z z=Z(a)/(Z(a)-b),l;acb_log1p(l.v,b.v,precision);
    return clog(Z(1)-Z(a)/b)*l-dilog(z*(Z(1)+b))+dilog(z);
}
Real inverse_triangle_full(const Real&a,const Real&b){
    pos(a,"positive triangle a");pos(b,"positive triangle b");
    Real alpha=sq(a-b)-2*(a+b)+1,beta=2*(a+b-1),root=rt(Real(5));pos(-alpha,"positive triangle Gram determinant");
    Real kappa=rt(-alpha),z0=beta/(2*root),q1=rt(sq(a-b)+4);
    Real delta=(4*sq(root)*alpha-sq(beta))/(2*root*(2*root*q1+2*sq(root)+beta));pos(-delta,"decreasing rational triangle chart");
    Z pole;arb_set(acb_imagref(pole.v),kappa.v);Z kernel=Z(delta)/(Z(z0)-pole),lk;acb_log1p(lk.v,kernel.v,precision);
    Z result=Z(2*real_log((root+1)/(root-1)))*lk;
    std::array<Real,2>linear={a-b-1,-a+b-1},square={b,a};
    for(int j=0;j<2;++j)for(long sign:{1L,-1L}){
        Real leading=Real(sign)-root,middle=beta-2*root*Real(sign)*linear[j],discriminant=4*rt(square[j]*(square[j]+4));
        for(long side:{1L,-1L}){Real r=(-middle+Real(side)*discriminant)/(2*leading),factor=delta/(z0-r);
            result=result+Z(sign)*pair_log_integral(factor,kernel);}
    }
    Real im;arb_set(im.v,acb_imagref(result.v));return -2*im/kappa;
}
Real triangle_log_full(const Real&a,const Real&b){
    Real cross=(a+b-1)/2,det=a*b-sq(cross);pos(det,"strict original triangle geometry");
    Real cx=b*(a-cross)/(2*det),cy=a*(b-cross)/(2*det),height=1+(a*cx+b*cy)/2;
    Real boundary=cx*edge_log_full(b)+cy*edge_log_full(a)+(1-cx-cy)*edge_log_full(Real(1));
    return boundary/2-Real(1)/2+height*inverse_triangle_full(a,b);
}
Real five_value(const Graph&g,bool positive,const std::array<Real,5>&p){
    Real s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],t=sq(s/(1-s)),A=1+t;
    Real den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den),co,si;
    arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real radial=rt(t*z),sqrt3=rt(Real(3));
    std::array<std::array<int,3>,3>parent={g.parent[0],g.parent[1],g.parent[0]};for(int j=0;j<3;++j)parent[2][j]-=g.parent[1][j];
    for(int j=0;j<3;++j)if(g.triangle[j][0]*g.triangle[j][0]!=(j<2?1:0))throw std::invalid_argument("proved loop partition required");
    const auto&fixed=g.triangle[2];if(fixed[1]*fixed[1]+fixed[2]*fixed[2]-fixed[1]*fixed[2]!=1)throw std::invalid_argument("unit fixed external invariant required");
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};Real V0;
    for(int j=0;j<3;++j)V0=V0+w[j]*Real(parent[j][0]*parent[j][0])*t;
    Real a0=bubble_A((1-x+V0)/x),K0=-(a0-bubble_A((1-x)/x)),T0=-edge_log_full(t)/2,total;
    for(long sign:{-1L,1L}){
        Real qx=Real(sign)*radial*co,qy=Real(sign)*radial*si,q1=qx,q2=(-qx+sqrt3*qy)/2;
        auto norm=[&](const std::array<int,3>&r){
            if(positive)return sq(Real(r[0])*qx+Real(r[1])-Real(r[2])/2)+sq(Real(r[0])*qy+sqrt3*Real(r[2])/2)+Real(r[0]*r[0])*t*(1-z);
            return Real(r[0]*r[0])*t+2*Real(r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real V,lr;for(int j=0;j<3;++j){V=V+w[j]*norm(parent[j]);Real host=1+norm(g.host[j]);pos(host,"positive host propagator");lr=lr-lp((host-A)/A);}
        Real DK=-(bubble_A((1-x+V)/x)-a0),DT=-triangle_log_full(norm(g.triangle[0]),norm(g.triangle[1]))-T0,ratio=real_exp(lr);
        total=total+(em1(lr)*K0*T0+ratio*(DK*T0+K0*DT+DK*DT))/2;
    }
    Real answer=2/(s*(1-s))*power3(t/A)*(1-z)*(1-x)*jac*total;
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite complete five-dimensional forest required");return answer;
}
Real five_controlled(const Graph&g,bool positive,const std::array<Real,5>&p){
    Real answer=five_value(g,positive,p),s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],t=sq(s/(1-s)),A=1+t;
    Real den=sq(x0)+sq(1-x0),x=sq(x0)/den,jac=2*x0*(1-x0)/sq(den),co,si;
    arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real rad=rt(t*z);std::array<Real,2>q={rad*co,rad*(-co+rt(Real(3))*si)/2};
    std::array<std::array<int,3>,3>parent={g.parent[0],g.parent[1],g.parent[0]};for(int j=0;j<3;++j)parent[2][j]-=g.parent[1][j];
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    Real V0,Cp,hsquare,hconstant;std::array<Real,2>pv,hv;
    for(int j=0;j<3;++j){const auto&a=parent[j];const auto&h=g.host[j];V0=V0+w[j]*Real(a[0]*a[0])*t;
        Cp=Cp+w[j]*Real(a[1]*a[1]+a[2]*a[2]-a[1]*a[2]);
        std::array<Real,2>one={Real(2*h[0]*h[1]),Real(2*h[0]*h[2])};hsquare=hsquare+norm_gram(one);
        hconstant=hconstant+Real(h[1]*h[1]+h[2]*h[2]-h[1]*h[2]);
        for(int k=0;k<2;++k){pv[k]=pv[k]+w[j]*Real(2*a[0]*a[k+1]);hv[k]=hv[k]+one[k];}}
    std::array<Real,2>a,b,tv,rv={-hv[0]/A,-hv[1]/A};
    for(int k=0;k<2;++k){a[k]=Real(g.triangle[0][0]*g.triangle[0][k+1]);b[k]=Real(g.triangle[1][0]*g.triangle[1][k+1]);tv[k]=a[k]+b[k];}
    Real h;arb_asinh(h.v,(rt(t)/2).v,precision);Real J0=4*h/rt(t*(t+4)),J2=2*(1+J0)/(t+4);
    Real U1=-(1-J0)/(2*t),Cpair=U1/2,Cfixed=-(J0-3*(1-J0)/t)/12,I2=(1-2*J0+J2)/sq(t),Ct;
    for(int j=0;j<3;++j){const auto&r=g.triangle[j];Ct=Ct+(j<2?Cpair:Cfixed)*Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);}
    Real Lt=q[0]*tv[0]+q[1]*tv[1],qa=q[0]*a[0]+q[1]*a[1],qb=q[0]*b[0]+q[1]*b[1];
    Real Lp=q[0]*pv[0]+q[1]*pv[1],R1=q[0]*rv[0]+q[1]*rv[1],hpoint;
    for(const auto&hh:g.host){Real one=2*Real(hh[0])*(Real(hh[1])*q[0]+Real(hh[2])*q[1]);hpoint=hpoint+sq(one);}
    Real R2=(sq(q[0]*hv[0]+q[1]*hv[1])+hpoint)/(2*sq(A))-hconstant/A;
    Real meanR2=t*(norm_gram(hv)+hsquare)/(12*sq(A))-hconstant/A;
    Real arg=(1-x+V0)/x,K0=-(bubble_A(arg)-bubble_A((1-x)/x)),B=bubble_value(arg);
    Real Kp=-x*B/sq(1+V0),Kpp=2*x*B/power3(1+V0)-bubble_derivative(arg)/sq(1+V0),T0=-edge_log_full(t)/2;
    Real point=R2*K0*T0+R1*(Kp*Lp*T0+K0*U1*Lt)+(Kp*Cp+Kpp*sq(Lp)/2)*T0+Kp*U1*Lp*Lt+K0*(Ct+I2*(sq(qa)+qa*qb+sq(qb))/3);
    Real mean=meanR2*K0*T0+t/6*(dot_gram(rv,pv)*Kp*T0+dot_gram(rv,tv)*K0*U1)+(Kp*Cp+Kpp*t*norm_gram(pv)/12)*T0+Kp*U1*t*dot_gram(pv,tv)/6+K0*(Ct+I2*t*(norm_gram(a)+dot_gram(a,b)+norm_gram(b))/18);
    Real common=2/(s*(1-s))*power3(t/A)*(1-x)*jac;
    answer=answer+common*(mean/2-(1-z)*point);
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite complete marginal control required");return answer;
}
Real five_mapped(const Graph&g,bool positive,bool matched,const std::array<Real,5>&p){
    if(!matched)return five_controlled(g,positive,p);
    Real x=sq(p[3])/(sq(p[3])+sq(1-p[3])),y=p[4],s=p[0];
    std::array<Real,3>w={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    std::array<int,3>loops={g.parent[0][0],g.parent[1][0],g.parent[0][0]-g.parent[1][0]};Real ell;
    for(int j=0;j<3;++j)ell=ell+w[j]*Real(loops[j]*loops[j]);pos(ell,"positive parent loop variance");
    Real h=1/rt(ell),den=1-s+h*s;auto pp=p;pp[0]=h*s/den;
    return h/sq(den)*five_controlled(g,positive,pp);
}
