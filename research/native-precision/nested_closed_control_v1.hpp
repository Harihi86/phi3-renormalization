// Included inside namespace nested_closed after the base full density.
Real dot_gram(const std::array<Real,2>&a,const std::array<Real,2>&b){return (a[0]-a[1]/2)*(b[0]-b[1]/2)+3*a[1]*b[1]/4;}
Real norm_gram(const std::array<Real,2>&a){return dot_gram(a,a);}
Real bubble_value(const Real&t){Real z=rt(t/(t+4)),h;arb_atanh(h.v,z.v,precision);return (t+4)/3*(h/z-1)-t/9;}
Real bubble_derivative(const Real&t){Real z=rt(t/(t+4)),h;arb_atanh(h.v,z.v,precision);return ((t-2)*h/z+2)/(3*t)-Real(5)/18;}
Real controlled_value(const Graph&g,bool positive,bool control,const std::array<Real,6>&p){
    Real answer=value(g,positive,p);if(!control)return answer;
    for(int j=0;j<3;++j)if(g.triangle[j][0]*g.triangle[j][0]!=(j<2?1:0))throw std::invalid_argument("proved triangle loop partition required for control");
    Real s=p[0],z=p[1],theta=p[2],x0=p[3],y=p[4],u=p[5],t=sq(s/(1-s)),A=1+t;
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
    for(int k=0;k<2;++k){a[k]=Real(g.triangle[0][0]*g.triangle[0][k+1]);b[k]=Real(g.triangle[1][0]*g.triangle[1][k+1]);tv[k]=u*(1-u)*(a[k]+b[k]);}
    Real Ct;for(int j=0;j<3;++j){const auto&r=g.triangle[j];Real weight=j<2?u*(1-u)/2:sq(1-u)/6;Ct=Ct+weight*Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);}
    Real Lt=q[0]*tv[0]+q[1]*tv[1],qa=q[0]*a[0]+q[1]*a[1],qb=q[0]*b[0]+q[1]*b[1];
    Real Lt2=4*sq(u*(1-u))*(sq(qa)+qa*qb+sq(qb))/3;
    Real tvnorm=4*sq(u*(1-u))*(norm_gram(a)+dot_gram(a,b)+norm_gram(b))/3;
    Real Lp=q[0]*pv[0]+q[1]*pv[1],R1=q[0]*rv[0]+q[1]*rv[1],hpoint;
    for(const auto&h:g.host){Real one=2*Real(h[0])*(Real(h[1])*q[0]+Real(h[2])*q[1]);hpoint=hpoint+sq(one);}
    Real R2=(sq(q[0]*hv[0]+q[1]*hv[1])+hpoint)/(2*sq(A))-hconstant/A;
    Real meanR2=t*(norm_gram(hv)+hsquare)/(12*sq(A))-hconstant/A;
    Real arg=(1-x+V0)/x,K0=-(bubble_A(arg)-bubble_A((1-x)/x)),B=bubble_value(arg);
    Real Kp=-x*B/sq(1+V0),Kpp=2*x*B/power3(1+V0)-bubble_derivative(arg)/sq(1+V0);
    Real VT0=u*(1-u)*t,T0=-lp(VT0),Tp=-Real(1)/(1+VT0),Tpp=Real(1)/sq(1+VT0);
    Real point=R2*K0*T0+R1*(Kp*Lp*T0+K0*Tp*Lt)+(Kp*Cp+Kpp*sq(Lp)/2)*T0+Kp*Tp*Lp*Lt+K0*(Tp*Ct+Tpp*Lt2/2);
    Real mean=meanR2*K0*T0+t/6*(dot_gram(rv,pv)*Kp*T0+dot_gram(rv,tv)*K0*Tp)+(Kp*Cp+Kpp*t*norm_gram(pv)/12)*T0+Kp*Tp*t*dot_gram(pv,tv)/6+K0*(Tp*Ct+Tpp*t*tvnorm/12);
    Real common=2/(s*(1-s))*power3(t/A)*(1-x)*jac*(1-u);
    answer=answer+common*(mean/2-(1-z)*point);
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite integrated angular control required");return answer;
}
