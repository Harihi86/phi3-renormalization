// Frozen scalar expressions copied verbatim from tested prefix kernel.
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::Z;
using fourier_e2_root::readq;using fourier_e2_root::interval;
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real lp(const Real&a){Real r;arb_log1p(r.v,a.v,precision);return r;}
Real real_log(const Real&a){Real r;arb_log(r.v,a.v,precision);return r;}
Real real_exp(const Real&a){Real r;arb_exp(r.v,a.v,precision);return r;}
Real em1(const Real&a){Real r;arb_expm1(r.v,a.v,precision);return r;}
Real power3(const Real&a){return a*sq(a);}
void pos(const Real&a,const char*msg){if(!arb_is_positive(a.v))throw std::runtime_error(msg);}
Z clog(const Z&a){Z r;acb_log(r.v,a.v,precision);return r;}
Z dilog(const Z&a){Z r;acb_polylog_si(r.v,2,a.v,precision);return r;}
Real atan_log(const Real&q){
    Z a;arb_sqrt_ui(acb_imagref(a.v),3,precision);Z sum;
    for(long b:{-1L,1L}){Z zz=-a/(Z(b)-a);Real u=1+Real(b)*q;pos(u,"positive marginal logarithm argument");
        sum=sum+(clog(Z(1)-a/Z(b))*Z(real_log(u))-dilog(zz*Z(u))+dilog(zz))/Z(b);}
    Real r;arb_set(r.v,acb_imagref(sum.v));return r/2;
}
Real bubble_A(const Real&t){
    if(arb_is_zero(t.v))return Real(0);pos(t,"nonnegative bubble primitive argument");
    Real q=rt(t/(t+4)),h,angle;pos(1-q,"bubble primitive finite q");
    arb_asinh(h.v,(rt(t)/2).v,precision);arb_atan(angle.v,(rt(Real(3))*q).v,precision);
    Real value=sq(h)/3+2*rt(Real(3))*(h*angle-atan_log(q))+4*q*h/(1+3*sq(q))-Real(17)/18*lp(t)-Real(32)/9*sq(q)/(1+3*sq(q));
    if(!arb_is_finite(value.v))throw std::runtime_error("finite bubble primitive required");return value;
}

Real weighted_C(const Real&t){
    if(arb_is_zero(t.v))return Real(0);pos(t,"nonnegative weighted primitive argument");
    Real h;arb_asinh(h.v,(rt(t)/2).v,precision);
    return h*rt(t*(t+4))*(t-2)/(3*(t+1))+4*sq(h)/3-Real(11)*t/18+lp(t)/18+8*t/(9*(t+1));
}
