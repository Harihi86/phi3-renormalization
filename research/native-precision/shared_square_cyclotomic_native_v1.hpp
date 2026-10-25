// Real six-periodic Li1/Li2/Li3 sum, with an analytic tail added as Arb error.
#include "shared_square_cyclotomic_coefficients.hpp"
Real cyclic_tail(const Real&w,const Real&logw,long count){
 Real n(count+1),l;arb_abs(l.v,logw.v);Real inv=1/n,inv2=sq(inv),inv3=inv*inv2,power;
 arb_pow_ui(power.v,w.v,count+1,precision);
 return power/(1-w)*(square_q(SHARED_CYCLIC_BOUND[0])*inv+
     square_q(SHARED_CYCLIC_BOUND[1])*(l*inv+inv2)+
     square_q(SHARED_CYCLIC_BOUND[2])*(sq(l)*inv+2*l*inv2+2*inv3));
}
struct CyclicResult{Real value,tail;long terms;bool series;};
CyclicResult shared_square_cyclic(const Real&w){
 pos(w,"positive cyclotomic argument required");pos(1-w,"cyclotomic argument below one required");
 if(!arb_le(w.v,(Real(1)/4).v))return {shared_square_P(w,false),Real(0),0,false};
 Real logw=real_log(w),target;arb_one(target.v);arb_mul_2exp_si(target.v,target.v,-precision-64);
 double estimate=arf_get_d(arb_midref(w.v),ARF_RND_UP);
 long count=1;
 if(estimate>0&&std::isfinite(estimate))count=std::max(1L,long(std::ceil((precision+96)*std::log(2.0)/(-std::log(estimate)))));
 if(count>8192)throw std::runtime_error("bounded cyclotomic term estimate required");
 Real error=cyclic_tail(w,logw,count);
 while(!arb_le(error.v,target.v)){
  if(count>=8192)throw std::runtime_error("cyclotomic tail did not reach declared absolute bound");
  count=std::min(count+8,8192L);error=cyclic_tail(w,logw,count);
 }
 static long cached_bits=0;
 static std::vector<std::array<Real,3>> coefficients;
 if(cached_bits!=precision){coefficients.clear();cached_bits=precision;}
 while(long(coefficients.size())<count){
  long n=long(coefficients.size())+1;Real inv=1/Real(n),inv2=sq(inv),inv3=inv*inv2;
  Real a=square_q(SHARED_CYCLIC[0][n%6]),b=square_q(SHARED_CYCLIC[1][n%6]),c=square_q(SHARED_CYCLIC[2][n%6]);
  coefficients.push_back({-a*inv+b*inv2-2*c*inv3,-b*inv+2*c*inv2,-c*inv});
 }
 std::array<Real,3> sums;
 for(long n=count;n>0;--n)for(long k=0;k<3;++k)sums[k]=sums[k]*w+coefficients[n-1][k];
 Real l2=sq(logw),value=w*(sums[0]+logw*sums[1]+l2*sums[2]);
 std::array<Real,3>powers={Real(1),logw,l2};
 for(long k=0;k<3;++k){const auto&part=SHARED_SQUARE_A[k];
  value=value+square_poly(part.numerator,w)/square_poly(part.denominator,w)*powers[k];
  value=value+square_q(part.real[0])*powers[k]*logw/Real(k+1);
 }
 arb_add_error(value.v,error.v);
 if(!arb_is_finite(value.v))throw std::runtime_error("finite bounded cyclotomic primitive required");
 return {value,error,count,true};
}
Real shared_square_P_fast(const Real&w){return shared_square_cyclic(w).value;}
