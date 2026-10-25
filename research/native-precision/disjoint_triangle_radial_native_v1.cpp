unsigned radial_series_terms(const Real&upper){
    if(arb_is_zero(upper.v))return 24;
    slong exponent=arf_abs_bound_lt_2exp_si(arb_midref(upper.v)),decay=std::max(slong(1),-exponent);
    return std::max(24u,unsigned((precision+decay-1)/decay)+8);
}
Real radial_phi(const Real&t){
    pos(t,"positive radial phi argument");Real pi;arb_const_pi(pi.v,precision);Real constant=sq(pi)/6-Real(5)/4;
    if(arb_is_one(t.v))return constant;
    Real delta=t-1,upper=arb_cw::abs_upper(delta),half=Real(1)/2;
    if(arb_le(upper.v,half.v)){
        unsigned count=radial_series_terms(upper);Real result=constant,power=delta;
        for(unsigned j=1;j<=count;++j){result=result+Real(j%2?2:-2)*power/(Real(long(j))*Real(long(j))*Real(long(j+1))*Real(long(j+2)));power=power*delta;}
        Real tail=2*Pow::power(upper,count+1)/(Real(long(count+1))*Real(long(count+1))*Real(long(count+2))*Real(long(count+3))*(1-upper));
        arb_cw::add_error(result,arb_cw::abs_upper(tail));return result;
    }
    Real dilog,arg=1-t;arb_polylog_si(dilog.v,2,arg.v,precision);
    return sq(pi)/6-dilog+t*(2-3*t)*logarithm(t)/(2*sq(1-t))-t/(2*(1-t));
}
std::array<Real,4> radial_moments(const Real&h){
    Real delta=h-1,upper=arb_cw::abs_upper(delta),half=Real(1)/2;std::array<Real,4>out;
    if(arb_le(upper.v,half.v)){
        unsigned count=radial_series_terms(upper);
        for(unsigned n=1;n<=4;++n){Real power(1),sum;for(unsigned j=0;j<=count;++j){sum=sum+power/Real(long(n+j));power=-power*delta;}
            Real tail=Pow::power(upper,count+1)/(Real(long(n+count+1))*(1-upper));arb_cw::add_error(sum,arb_cw::abs_upper(tail));out[n-1]=sum;}
    }else{out[0]=logarithm(h)/delta;for(int n=2;n<=4;++n)out[n-1]=(Real(1)/Real(n-1)-out[n-2])/delta;}
    return out;
}
Real radial_primitive(const Real&q,const Real&h,const Real&de,const Real&ge,const Real&mean_de,bool control){
    pos(h,"positive nonempty child variance");Real a=1+q,b=1+h*a+de,c=h+ge,variance=h*de-ge;
    // The exact Schur identity proves variance=|c_vector-h*b_vector|^2.
    // Intersecting its interval with the nonnegative axis is therefore valid.
    if(arb_is_negative(variance.v))throw std::runtime_error("inconsistent exact Schur variance");
    arb_nonnegative_part(variance.v,variance.v);
    Real disc=sq(1-h*a+de)+4*a*variance,s=root(disc),plus=(b+s)/(2*a),minus=2*c/(b+s),la=stable_log1p(q);
    Real result=la*logarithm(h)+sq(la)/2+3*la/2-radial_phi(plus/h)-radial_phi(minus/h)+radial_phi(1/h)+radial_phi(Real(1));
    if(control){auto moments=radial_moments(h);result=result+(Real(1)/4-q)*h*(moments[0]+moments[1]+moments[2])-(mean_de-de)*moments[3];}
    return result;
}
std::array<Real,3> common_radial_coefficients(const Graph&g,const std::vector<std::pair<int,int>>&rows,const std::vector<Real>&alpha,const std::vector<std::array<Real,3>>&child,const std::array<Real,3>&y,int active,bool mean){
    std::array<Real,3>out;static const int ij[3][2]={{0,1},{0,2},{1,2}};
    for(const auto&entry:rows){int mask=entry.first,degree=0;Real monomial(entry.second);bool use=true;
        if(mean){int d=pop(mask&7),factorial=1;for(int j=2;j<=d+2;++j)factorial*=j;monomial=monomial*2/factorial;}
        else for(int e=0;e<3;++e)if(mask&(1<<e))monomial=monomial*alpha[e];
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(mask&(1<<(3+3*i+j))){
            if(!(active&(1<<i))){use=false;break;}
            ++degree;monomial=monomial*y[i]*child[i][ij[j][0]]*child[i][ij[j][1]];}
        if(use){if(degree>2)throw std::runtime_error("unexpected radial polynomial degree");out[degree]=out[degree]+monomial;}
    }
    return out;
}
Real gram_norm(const std::array<Real,2>&v){return sq(v[0]-v[1]/2)+3*sq(v[1])/4;}
Real gram_dot(const std::array<Real,2>&v,const std::array<Real,2>&w){return (v[0]-v[1]/2)*(w[0]-w[1]/2)+3*v[1]*w[1]/4;}
Real radial_density(const Graph&g,bool original,bool control,const std::array<Rat,11>&input){
    if(g.id!="L4E3G521"||g.L!=1||g.k!=3)throw std::invalid_argument("integrated radial mode requires the proved three-child graph");
    std::array<Real,11>p;for(int j=0;j<11;++j)p[j]=real(input[j]);auto hostx=simplex(p[0],p[1]),y=simplex(p[8],p[9]);
    std::vector<Real>alpha={hostx[0],hostx[1],hostx[2]};std::vector<std::array<Real,3>>child;
    for(int i=0;i<3;++i)child.push_back(simplex(p[2+2*i],p[3+2*i]));
    Real host[3][3];for(int e=0;e<3;++e)for(int i=0;i<3;++i)for(int j=0;j<3;++j)host[i][j]=host[i][j]+alpha[e]*Real(g.routes[e][i]*g.routes[e][j]);
    std::array<Real,2>hb={host[0][1],host[0][2]},center;
    for(int j=0;j<2;++j)for(int e=0;e<3;++e)center[j]=center[j]+Real(g.routes[e][0]*g.routes[e][j+1])/3;
    Real hostq=host[1][1]+host[2][2]-(host[1][2]+host[2][1])/2-gram_norm(hb),total;
    static const int ij[3][2]={{0,1},{0,2},{1,2}};
    for(int active=1;active<8;++active){Real q,h,de,ge,mean_de;
        if(original){Real variance[3][3];
            for(int c=0;c<3;++c)if(active&(1<<c))for(int edge=0;edge<3;++edge)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                variance[i][j]=variance[i][j]+y[c]*child[c][ij[edge][0]]*child[c][ij[edge][1]]*Real(g.delta[c][edge][i]*g.delta[c][edge][j]);
            h=variance[0][0];std::array<Real,2>cross={variance[0][1],variance[0][2]};Real energy=variance[1][1]+variance[2][2]-(variance[1][2]+variance[2][1])/2;
            q=hostq;de=energy-2*gram_dot(hb,cross)+h*gram_norm(hb);ge=h*energy-gram_norm(cross);
            mean_de=energy-2*gram_dot(center,cross)+h*gram_norm(center)+h/12;
        }else{auto u=common_radial_coefficients(g,g.U,alpha,child,y,active,false),f=common_radial_coefficients(g,g.F,alpha,child,y,active,false);
            h=u[1];q=f[0];de=f[1]-h*q;ge=f[2];if(control){auto mean=common_radial_coefficients(g,g.F,alpha,child,y,active,true);mean_de=mean[1]-h/4;}}
        total=total+Real((3-pop(active))%2?-1:1)*radial_primitive(q,h,de,ge,mean_de,control);
    }
    total=total/(32*y[0]*y[1]*y[2]);if(!arb_is_finite(total.v))throw std::runtime_error("finite complete analytic radial density required");return total;
}
