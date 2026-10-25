// One bubble, three host loops, two physical external momenta, nine complete sectors.
#include "exact_rational_sector_native_v1.hpp"
#include "single_weighted_e3_tables.hpp"
Rat single_polynomial(const std::vector<std::array<int,2>>&terms,const std::array<Rat,9>&a){
 Rat total;for(const auto&term:terms){Rat value{long(term[1])};for(long e=0;e<9;++e)if(term[0]&(1L<<e))value=value*a[e];total=total+value;}return total;
}
std::array<Rat,2> single_exact_uf(const State&s,const std::array<Rat,9>&a,bool positive){
 if(positive)return {single_polynomial(s.U,a),single_polynomial(s.F,a)};
 std::array<std::array<Rat,5>,5>m;
 for(long e=0;e<9;++e)for(long i=0;i<5;++i)for(long j=0;j<5;++j)m[i][j]=m[i][j]+a[e]*Rat(long(s.routes[e][i]*s.routes[e][j]));
 std::array<std::array<Rat,3>,3>loop;for(long i=0;i<3;++i)for(long j=0;j<3;++j)loop[i][j]=m[i][j];
 auto minor=[&](long row,long col){std::array<std::array<Rat,4>,4>b;std::array<long,4>ri={0,1,2,3+row},ci={0,1,2,3+col};
  for(long i=0;i<4;++i)for(long j=0;j<4;++j)b[i][j]=m[ri[i]][ci[j]];return rational_determinant(b);};
 return {rational_determinant(loop),minor(0,0)+minor(1,1)-(minor(0,1)+minor(1,0))/Rat(2)};
}
struct SingleTerm{Rat t1,coefficient;};
struct SingleSector{Rat t0,weight,constant;std::vector<SingleTerm>terms;};
std::vector<SingleSector> prepare_single(const Graph&g,const std::array<Rat,9>&p,bool positive,bool fourth){
 const auto&t=single_sector_table(g.id,fourth);std::vector<SingleSector>sectors;
 long weighted=g.states.at(0).weighted;
 for(long orbit=0;orbit<9;++orbit){Rat selector=p[0]+Rat(orbit)/Rat(9);if(!(selector<Rat(1)))selector=selector-Rat(1);
  auto mapped=rational_sector_point(t,selector,p);Rat mass;for(const auto&a:mapped.alpha)mass=mass+a;Rat d=mapped.alpha[weighted];
  SingleSector sector;sector.t0=(mass-d)/d;sector.weight=mapped.weight/Rat(9);
  for(const auto&state:g.states){
   if(long(state.weighted)!=weighted)throw std::runtime_error("fixed original bubble-weighted edge required");
   auto pol=single_exact_uf(state,mapped.alpha,positive);
   if(!(Rat(0)<pol[0])||!(Rat(0)<pol[1]))throw std::runtime_error("positive exact U and F required");
   Rat coefficient=Rat(long(state.sign))/rat_power(pol[0],3),t1=sector.t0+pol[1]/(pol[0]*d);
   sector.constant=sector.constant+coefficient;sector.terms.push_back({t1,coefficient});
  }
  sectors.push_back(sector);
 }
 return sectors;
}
Real single_value(const Graph&g,const std::array<Rat,9>&p,bool positive,bool fourth,long requested_bits){
 auto sectors=prepare_single(g,p,positive,fourth);std::string last_error;
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{Real total;
   for(const auto&sector:sectors){Real subtotal;
    if(!fmpq_is_zero(sector.constant.v))subtotal=to_real(sector.constant)*bubble_A(to_real(sector.t0));
    for(const auto&term:sector.terms)subtotal=subtotal-to_real(term.coefficient)*bubble_A(to_real(term.t1));
    total=total+to_real(sector.weight)*subtotal;
   }
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){std::cerr<<"single-sector-working-bits "<<working<<'\n';return total;}
   last_error="complete single-bubble E3 forest missed80relativebits";
  }catch(const std::runtime_error&e){last_error=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("bounded single-bubble precision exhausted: "+last_error);
}
