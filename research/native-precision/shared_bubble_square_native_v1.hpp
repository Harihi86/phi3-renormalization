// Exact rational coefficients come from the actual Hermite derivative proof.
using SquareQ=std::array<long,2>;
struct SquarePart{std::vector<SquareQ>numerator,denominator;std::array<SquareQ,3>real;std::array<SquareQ,2>rho;};
#include "shared_bubble_square_coefficients.hpp"
Real square_q(const SquareQ&v){return Real(v[0])/Real(v[1]);}
Real square_poly(const std::vector<SquareQ>&coefficients,const Real&w){Real out;for(const auto&c:coefficients)out=out*w+square_q(c);return out;}
Real square_real(const Z&z){Real out;arb_set(out.v,acb_realref(z.v));return out;}
Z square_polylog(long order,const Z&z){Z out;acb_polylog_si(out.v,order,z.v,precision);return out;}
std::array<Z,3> square_phi(const Z&pole,const Real&w,const Real&logw){
 Z z=Z(w)/pole,ell=clog(Z(1)-z),li2=square_polylog(2,z),li3=square_polylog(3,z),l(logw);
 return {ell,l*ell+li2,l*l*ell+Z(2)*l*li2-Z(2)*li3};
}
Real shared_square_P(const Real&w,bool at_one){
 Real logw=at_one?Real(0):real_log(w),l2=sq(logw),total;
 std::array<Real,3>powers={Real(1),logw,l2};
 for(long k=0;k<3;++k){const auto&s=SHARED_SQUARE_A[k];
  total=total+square_poly(s.numerator,w)/square_poly(s.denominator,w)*powers[k];
  total=total+square_q(s.real[0])*powers[k]*logw/Real(k+1);
 }
 for(long index=1;index<3;++index){
  bool nonzero=false;for(long k=0;k<3;++k)nonzero=nonzero||SHARED_SQUARE_A[k].real[index][0]!=0;
  if(!nonzero)continue;
  if(index==1&&at_one){
   if(SHARED_SQUARE_A[0].real[1][0]||SHARED_SQUARE_A[1].real[1][0])throw std::runtime_error("finite exact normalization required");
   Real zeta;arb_zeta_ui(zeta.v,3,precision);total=total-2*square_q(SHARED_SQUARE_A[2].real[1])*zeta;
  }else{
   auto phi=square_phi(Z(index==1?1L:-1L),w,logw);
   for(long k=0;k<3;++k)total=total+square_q(SHARED_SQUARE_A[k].real[index])*square_real(phi[k]);
  }
 }
 Z rho;arb_set_si(acb_realref(rho.v),1);arb_sqrt_ui(acb_imagref(rho.v),3,precision);rho=rho/Z(2);
 auto phi=square_phi(rho,w,logw);
 for(long k=0;k<3;++k){const auto&s=SHARED_SQUARE_A[k];Z coefficient;
  arb_set(acb_realref(coefficient.v),square_q(s.rho[0]).v);
  Real imaginary=rt(Real(3))*square_q(s.rho[1]);arb_set(acb_imagref(coefficient.v),imaginary.v);
  total=total+2*square_real(coefficient*phi[k]);
 }
 return total;
}
Real shared_square_A(const Real&t){
 if(arb_is_zero(t.v))return Real(0);pos(t,"positive shared-square primitive argument");
 static long cached_bits=0;static Real normalization;
 if(cached_bits!=precision){normalization=shared_square_P(Real(1),true);cached_bits=precision;}
 Real w=4/sq(rt(t+4)+rt(t));pos(w,"positive rationalized shared-square coordinate");pos(1-w,"shared-square coordinate below one");
 Real out=shared_square_P(w,false)-normalization;
 if(!arb_is_finite(out.v))throw std::runtime_error("finite shared-square primitive required");return out;
}
