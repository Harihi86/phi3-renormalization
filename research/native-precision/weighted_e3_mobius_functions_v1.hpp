// Complete E3 forest, both external momenta, unchanged shared/distinct spectral kernels.
#include "exact_rational_sector_native_v1.hpp"
#include "weighted_e3_mobius_tables_v1.hpp"
Rat e3_rational_polynomial(const std::vector<std::array<int,2>>&terms,const std::array<Rat,6>&a){
 Rat sum;for(const auto&term:terms){Rat product{long(term[1])};for(long e=0;e<6;++e)if(term[0]&(1L<<e))product=product*a[e];sum=sum+product;}return sum;
}
std::array<Rat,2> e3_exact_uf(const State&s,const std::array<Rat,6>&a,bool positive){
 if(positive)return {e3_rational_polynomial(s.U,a),e3_rational_polynomial(s.F,a)};
 std::array<std::array<Rat,4>,4>m;
 for(long e=0;e<6;++e)for(long i=0;i<4;++i)for(long j=0;j<4;++j)m[i][j]=m[i][j]+a[e]*Rat(long(s.routes[e][i]*s.routes[e][j]));
 auto minor=[&](long row,long col){std::array<std::array<Rat,3>,3>b;std::array<long,3>ri={0,1,2+row},ci={0,1,2+col};
  for(long i=0;i<3;++i)for(long j=0;j<3;++j)b[i][j]=m[ri[i]][ci[j]];return rational_determinant(b);};
 Rat U=m[0][0]*m[1][1]-m[0][1]*m[1][0];
 return {U,minor(0,0)+minor(1,1)-(minor(0,1)+minor(1,0))/Rat(2)};
}
struct E3SectorTerm{Rat U,q;long sign;};
struct E3Sector{Rat mass,d,e,weight;std::vector<E3SectorTerm>terms;};
std::vector<E3Sector> prepare_e3_sectors(const Graph&g,const std::array<Rat,7>&p,bool positive,bool fourth){
 const auto&t=e3_mobius_table(g.id,fourth);std::array<Rat,6>point;for(long j=0;j<6;++j)point[j]=p[j];std::vector<E3Sector>sectors;
 const auto&weights=g.states.at(0).weighted;
 if(weights.size()!=(g.shared?1UL:2UL))throw std::runtime_error("proved E3 weighted edges required");
 for(long orbit=0;orbit<6;++orbit){Rat selector=p[0]+Rat(orbit)/Rat(6);if(!(selector<Rat(1)))selector=selector-Rat(1);auto mapped=rational_sector_point(t,selector,point);
  E3Sector sector;for(const auto&a:mapped.alpha)sector.mass=sector.mass+a;sector.d=mapped.alpha[weights[0]];if(!g.shared)sector.e=mapped.alpha[weights[1]];sector.weight=mapped.weight/Rat(6);
  for(const auto&state:g.states){if(state.weighted!=weights)throw std::runtime_error("fixed original weighted edges required");auto uf=e3_exact_uf(state,mapped.alpha,positive);
   if(!(Rat(0)<uf[0])||!(Rat(0)<uf[1]))throw std::runtime_error("positive exact E3 U and F required");sector.terms.push_back({uf[0],uf[1]/uf[0],long(state.sign)});}
  sectors.push_back(sector);
 }
 return sectors;
}
Real e3_sector_value(const Graph&g,const std::array<Rat,7>&p,bool positive,bool fourth,long requested_bits){
 auto sectors=prepare_e3_sectors(g,p,positive,fourth);Rat v=p[6],r=Rat(1)-rat_power(v,2),tau=(Rat(1)-rat_power(r,2))/Rat(4);
 Rat spectral=g.shared?Rat(1):v*rat_power(r,4)/Rat(6);std::string last_error;
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{Real total;
   for(const auto&sector:sectors){Real subtotal;
    if(g.shared){
     for(const auto&term:sector.terms){Rat a=sector.mass+v*term.q;Real b=bubble_B(to_real((a-sector.d)/sector.d));
      Rat coefficient=Rat(-term.sign)*rat_power(sector.d,2)*term.q/(rat_power(a,3)*rat_power(term.U,3));subtotal=subtotal+to_real(coefficient)*sq(b);}
    }else{
     Real d=to_real(sector.d),e=to_real(sector.e);auto base=prepare(to_real(sector.mass),d),shifted=prepare(to_real(sector.mass+sector.e*(Rat(1)/tau-Rat(1))),d);
     Real c=to_real((tau-Rat(2))/rat_power(Rat(1)-tau,2)),h=to_real(Rat(1)/(Rat(1)-tau)),k=to_real(Rat(1)/(tau*rat_power(Rat(1)-tau,2)));
     for(const auto&term:sector.terms){auto j=J_mass(base,to_real(term.q)),other=J_mass(shifted,to_real(term.q));subtotal=subtotal+to_real(Rat(term.sign)/rat_power(term.U,3))*(c*j[0]-h*e*j[1]+k*other[0]);}
    }
    total=total+to_real(sector.weight*spectral)*subtotal;
   }
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){std::cerr<<"sector-working-bits "<<working<<'\n';return total;}
   last_error="complete E3 forest missed80relativebits";
  }catch(const std::runtime_error&e){last_error=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("bounded E3 sector precision exhausted: "+last_error);
}
