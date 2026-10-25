// The full shared-edge auxiliary integral is eliminated before sampling.
struct ClosedSharedTerm{Rat t1,coefficient;};
struct ClosedSharedSector{Rat t0,constant,weight;std::vector<ClosedSharedTerm>terms;};
std::vector<ClosedSharedSector> prepare_shared_closed(const Graph&g,const std::array<Rat,6>&p,bool positive,bool fourth){
 if(!g.shared||g.states.at(0).weighted.size()!=1)throw std::invalid_argument("proved shared-edge E3 graph required");
 const auto&t=e3_mobius_table(g.id,fourth);std::vector<ClosedSharedSector>sectors;const auto&weights=g.states.at(0).weighted;
 for(long orbit=0;orbit<6;++orbit){Rat selector=p[0]+Rat(orbit)/Rat(6);if(!(selector<Rat(1)))selector=selector-Rat(1);
  auto mapped=rational_sector_point(t,selector,p);Rat mass;for(const auto&a:mapped.alpha)mass=mass+a;
  Rat d=mapped.alpha[weights[0]];ClosedSharedSector sector;sector.t0=(mass-d)/d;sector.weight=mapped.weight/Rat(6);
  for(const auto&state:g.states){
   if(state.weighted!=weights)throw std::runtime_error("original shared edge identity required");auto uf=e3_exact_uf(state,mapped.alpha,positive);
   if(!(Rat(0)<uf[0])||!(Rat(0)<uf[1]))throw std::runtime_error("positive original E3 U and F required");
   Rat coefficient=Rat(long(state.sign))/rat_power(uf[0],3);
   sector.constant=sector.constant+coefficient;sector.terms.push_back({sector.t0+uf[1]/(uf[0]*d),coefficient});
  }
  sectors.push_back(sector);
 }
 return sectors;
}
Real shared_closed_value(const Graph&g,const std::array<Rat,6>&p,bool positive,bool fourth,long requested_bits){
 auto sectors=prepare_shared_closed(g,p,positive,fourth);std::string last_error;
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{Real total;
   for(const auto&sector:sectors){Real subtotal;
    if(!fmpq_is_zero(sector.constant.v))subtotal=to_real(sector.constant)*shared_square_A(to_real(sector.t0));
    for(const auto&term:sector.terms)subtotal=subtotal-to_real(term.coefficient)*shared_square_A(to_real(term.t1));
    total=total+to_real(sector.weight)*subtotal;
   }
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){std::cerr<<"shared-closed-working-bits "<<working<<'\n';return total;}
   last_error="complete shared E3 closed forest missed80relativebits";
  }catch(const std::runtime_error&e){last_error=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("bounded shared E3 precision exhausted: "+last_error);
}
