// Exact tails: all poles are nonreal; the upper endpoint uses a finite limit.
Z tail_power(Z a,int n){if(n<0)return Z(1)/tail_power(a,-n);Z out(1);while(n){if(n&1)out=out*a;a=a*a;n>>=1;}return out;}
Z tail_rational_moment(const Real&r,const Z&z,int j,int power){
    Z out;for(int l=0;l<=power;++l){long binomial=1;for(int k=1;k<=l;++k)binomial=binomial*(power-k+1)/k;
        int exponent=l-j+1;Z primitive=exponent==0?clog(Z(1)-z)-clog(Z(r)-z):(tail_power(Z(1)-z,exponent)-tail_power(Z(r)-z,exponent))/Z(exponent);
        out=out+Z(binomial)*tail_power(z,power-l)*primitive;}return out;
}
std::array<Real,3> bubble_tails(const Real&a){
    if(!arb_is_nonnegative(a.v))throw std::invalid_argument("nonnegative bubble tail argument");
    Real r=rt(a/(a+4));pos(1-r,"strict rationalized tail endpoint");Z z;arb_set(acb_imagref(z.v),(Real(1)/rt(Real(3))).v);
    std::array<Z,6>integrals,L,J,pieces;
    for(long sign:{-1L,1L}){
        Z A=Z(1)+Z(sign)*z,top=clog(A)*clog(Z(1)-z)-(sign==-1?dilog(Z(1)):dilog(-Z(sign)*(Z(1)-z)/A));
        Z low=clog(A)*clog(Z(r)-z)-dilog(-Z(sign)*(Z(r)-z)/A);
        integrals[1]=integrals[1]+Z(sign)*(top-low)/Z(2);
    }
    Real ar;arb_atanh(ar.v,r.v,precision);L[1]=clog(Z(1)-z)-clog(Z(r)-z);J[0]=Z(real_log(2/(1+r)));
    for(int k=1;k<=4;++k){
        if(k>1)L[k]=(tail_power(Z(1)-z,1-k)-tail_power(Z(r)-z,1-k))/Z(1-k);
        J[k]=(L[k]-J[k-1])/(Z(1)+z);
        Z vv=(tail_power(Z(r)-z,-k)-tail_power(Z(1)-z,-k))/Z(-k),sum;
        for(int m=1;m<=k;++m)sum=sum+tail_power(Z(1)-z,m-k-1)*J[m];
        integrals[k+1]=-vv*Z(ar)+sum/Z(k);
    }
    for(int j=1;j<=5;++j)pieces[j]=integrals[j]-tail_rational_moment(r,z,j,1)-tail_rational_moment(r,z,j,3)/Z(3);
    auto coefficient=[](long rn,long rd,long in,long id){Z c(Real(rn)/Real(rd));arb_set(acb_imagref(c.v),(rt(Real(3))*Real(in)/Real(id)).v);return c;};
    std::array<std::array<Z,5>,3>cs={{{coefficient(0,1,-2,3),coefficient(-2,3,0,1),coefficient(0,1,4,27),Z(0),Z(0)},
        {coefficient(0,1,-14,27),coefficient(-14,27,0,1),coefficient(0,1,4,27),coefficient(8,81,0,1),Z(0)},
        {coefficient(0,1,-4,9),coefficient(-4,9,0,1),coefficient(0,1,32,243),coefficient(8,81,0,1),coefficient(0,1,-16,729)}}};
    std::array<Real,3>answer;for(int n=0;n<3;++n){Z sum;for(int j=1;j<=n+3;++j)sum=sum+cs[n][j-1]*pieces[j];
        arb_mul_2exp_si(answer[n].v,acb_realref(sum.v),1);pos(answer[n],"positive exact bubble tail");}
    return answer;
}
std::array<Real,4> parent_integrated_moments(const Real&xi,const Real&t){
    pos(xi,"positive xi");pos(1-xi,"xi below one");Real X=1-xi,d=1+xi*X*t,a=d/X-1,b=1/X-1;
    auto tails=bubble_tails(a),base=bubble_tails(b);
    Real I0=-X*(bubble_A(a)-bubble_A(b))-d*tails[0]+base[0],I1=-d*tails[1],Ic=I1+sq(d)*tails[2],I2=d*bubble_value(a)/sq(sq(1+a))-2*d*tails[2];
    return {I0,I1,Ic,I2};
}
Real parent_integrated_control_mean(const TwoBubbleGraph&g,const Real&t,const Real&xi){
    if(g.parent[0][0]!=-1||g.parent[1][0]!=0||g.parent[0][1]!=g.parent[1][1]||g.parent[0][2]!=g.parent[1][2]
       ||g.parent[0][1]*g.parent[0][1]+g.parent[0][2]*g.parent[0][2]-g.parent[0][1]*g.parent[0][2]!=1)throw std::invalid_argument("proved graph-specific parent variance routing required");
    Real A=1+t,bt=bubble_value(t),J0=bt/A;std::array<Real,2>hv;Real hsquare,hconstant;
    for(const auto&h:g.host){std::array<Real,2>one={Real(2*h[0]*h[1]),Real(2*h[0]*h[2])};hsquare=hsquare+norm_gram(one);
        hconstant=hconstant+Real(h[1]*h[1]+h[2]*h[2]-h[1]*h[2]);for(int k=0;k<2;++k)hv[k]=hv[k]+one[k];}
    std::array<Real,2>rv={-hv[0]/A,-hv[1]/A},bv={Real(2*g.bubble[1]),Real(2*g.bubble[2])},pv={Real(-2*g.parent[0][1])*xi,Real(-2*g.parent[0][2])*xi};
    Real Cb=Real(g.bubble[1]*g.bubble[1]+g.bubble[2]*g.bubble[2]-g.bubble[1]*g.bubble[2]);
    Real meanR2=t*(norm_gram(hv)+hsquare)/(12*sq(A))-hconstant/A;
    Real hp;arb_asinh(hp.v,(rt(t)/2).v,precision);Real edge=4*hp/rt(t*(t+4)),bp=bubble_derivative(t),bpp=1/(6*t)-(1-edge)/sq(t);
    Real Jp=bp/A-bt/sq(A),Jpp=bpp/A-2*bp/sq(A)+2*bt/power3(A);
    auto moments=parent_integrated_moments(xi,t);Real I0=moments[0],I1=moments[1],Ic=moments[2],I2=moments[3];
    return meanR2*I0*J0+t/6*(dot_gram(rv,pv)*I1*J0+dot_gram(rv,bv)*I0*Jp)
        +(Ic+I2*t*norm_gram(pv)/12)*J0+I1*Jp*t*dot_gram(pv,bv)/6+I0*(Jp*Cb+Jpp*t*norm_gram(bv)/12);
}
