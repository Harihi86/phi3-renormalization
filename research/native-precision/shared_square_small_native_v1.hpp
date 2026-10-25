// A(t), normalized at zero, from the proved symmetric-square Taylor recurrence.
// Coefficient arithmetic and the analytic remainder both carry Arb enclosures.
#include "shared_square_small_operator.hpp"
struct SmallSquareResult{Real value,tail;long terms;bool series;long working_bits;};
Real small_square_lower(const Real&t){return t*sq(sq(t))/(18000*sq(1+t/4)*(1+t)*sq(1+t));}
Real small_square_tail(const Real&t,long count){
 Real power;arb_pow_ui(power.v,t.v,count+5,precision);
 return power*(Real(count+2)-Real(count+1)*t)/(4050*sq(1-t));
}
SmallSquareResult shared_square_small(const Real&t){
 if(arb_is_zero(t.v))return {Real(0),Real(0),0,true,precision};
 pos(t,"positive small-square argument required");
 if(!arb_le(t.v,(Real(1)/16).v))throw std::invalid_argument("proved small-square domain required");
 Real target=small_square_lower(t);arb_mul_2exp_si(target.v,target.v,-precision-32);
 double approximate=arf_get_d(arb_midref(t.v),ARF_RND_UP);long count=1;
 if(approximate>0&&std::isfinite(approximate))count=std::max(1L,long(std::ceil((precision+64)*std::log(2.0)/(-std::log(approximate)))));
 if(count>4096)throw std::runtime_error("bounded small-square term estimate required");
 Real tail=small_square_tail(t,count);
 while(!arb_le(tail.v,target.v)){
  if(count>=4096)throw std::runtime_error("small-square analytic tail bound exhausted");
  count=std::min(count+8,4096L);tail=small_square_tail(t,count);
 }
 static long cached_bits=0;static std::vector<Real> f;
 if(cached_bits!=precision){f.assign(5,Real(0));f[4]=Real(1)/3600;cached_bits=precision;}
 while(long(f.size())<count+4){
  long n=long(f.size()),leading=0;Real previous;
  for(long j=0;j<4;++j)for(const auto&row:SMALL_SQUARE_OPERATOR[j]){
   long index=n+j-row[0];if(index<j)continue;
   long factor=row[1];for(long k=0;k<j;++k)factor*=index-k;
   if(index==n)leading+=factor;
   else{if(index<0||index>=n)throw std::runtime_error("strict Taylor recurrence ordering required");previous=previous+Real(factor)*f.at(index);}
  }
  if(!leading)throw std::runtime_error("nonzero Taylor recurrence pivot required");
  f.push_back(-previous/Real(leading));
 }
 Real value;
 for(long n=count-1;n>=0;--n)value=value*t+f.at(n+4)/Real(n+5);
 value=value*t*sq(sq(t));arb_add_error(value.v,tail.v);
 pos(value,"positive normalized squared-bubble primitive required");
 if(!arb_is_finite(value.v))throw std::runtime_error("finite Taylor enclosure required");
 return {value,tail,count,true,precision};
}
