// Use the proved scalar kernel and exact U/F code; sum the conditional charts explicitly.
struct PrefixStratum{std::vector<long>prefix;Rat left,mass;};
std::vector<PrefixStratum> prefix_strata(const MobiusTable&t,long depth){
 if(depth!=1&&depth!=2)throw std::invalid_argument("one or two exact prefix choices required");
 std::vector<PrefixStratum>out;Rat left1;
 for(long first=0;first<8;++first){long remaining=255^(1L<<first);Rat mass1=Rat(t.weights[remaining])/Rat(t.normalization);
  if(depth==1)out.push_back({{first},left1,mass1});
  else{long total=0;for(long second=0;second<8;++second)if(remaining&(1L<<second))total+=t.weights[remaining^(1L<<second)];Rat left2;
   for(long second=0;second<8;++second)if(remaining&(1L<<second)){
    Rat mass2=Rat(t.weights[remaining^(1L<<second)])/Rat(total);
    out.push_back({{first,second},left1+mass1*left2,mass1*mass2});left2=left2+mass2;
   }
  }
  left1=left1+mass1;
 }
 return out;
}
std::vector<MobiusSector> prepare_prefix(const Graph&g,const std::array<Rat,8>&p,bool positive,long depth){
 const auto&t=mobius_table(g.id,false);std::vector<MobiusSector>sectors;
 for(const auto&stratum:prefix_strata(t,depth)){
  auto decoded=mobius_point(t,stratum.left+stratum.mass*p[0],p);
  for(long j=0;j<depth;++j)if(decoded.order[j]!=stratum.prefix[j])throw std::runtime_error("exact prefix identity required");
  Rat mass;for(const auto&a:decoded.alpha)mass=mass+a;long weighted=g.states.at(0).weighted;Rat d=decoded.alpha[weighted];
  MobiusSector sector;sector.t0=(mass-d)/d;sector.weight=decoded.weight*stratum.mass;
  for(const auto&state:g.states){
   if(long(state.weighted)!=weighted)throw std::runtime_error("original fixed weighted edge required");
   auto uf=exact_uf(state,decoded.alpha,positive);
   if(!(Rat(0)<uf[0])||!(Rat(0)<uf[1]))throw std::runtime_error("positive exact original U and F required");
   Rat delta=uf[1]/(uf[0]*d);sector.terms.push_back({sector.t0+delta,Rat(long(state.sign))*d/rat_power(uf[0],3)});
  }
  sectors.push_back(sector);
 }
 return sectors;
}
Real prefix_value(const Graph&g,const std::array<Rat,8>&p,bool positive,long depth,long requested_bits){
 auto sectors=prepare_prefix(g,p,positive,depth);std::string last_error;
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{Real total;
   for(const auto&sector:sectors){Real t0=to_real(sector.t0),a0=bubble_A(t0),c0=weighted_C(t0),subtotal;
    for(const auto&term:sector.terms){Real t1=to_real(term.t1);subtotal=subtotal+to_real(term.coefficient)*(t1*(bubble_A(t1)-a0)-(weighted_C(t1)-c0));}
    total=total+to_real(sector.weight)*subtotal;
   }
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){std::cerr<<"prefix-working-bits "<<working<<'\n';return total;}
   last_error="complete prefix forest missed80relativebits";
  }catch(const std::runtime_error&e){last_error=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("bounded prefix precision exhausted: "+last_error);
}
void prefix_map_output(const Graph&g,const std::array<Rat,8>&p,long depth){
 const auto&t=mobius_table(g.id,false);std::cout<<'[';long count=0;
 for(const auto&s:prefix_strata(t,depth)){auto d=mobius_point(t,s.left+s.mass*p[0],p);
  if(count++)std::cout<<',';std::cout<<"{\"prefix\":[";for(long j=0;j<depth;++j){if(j)std::cout<<',';std::cout<<s.prefix[j];}
  std::cout<<"],\"left\":\""<<s.left.str()<<"\",\"mass\":\""<<s.mass.str()<<"\",\"alpha\":[";
  for(long j=0;j<8;++j){if(j)std::cout<<',';std::cout<<'\"'<<d.alpha[j].str()<<'\"';}
  std::cout<<"],\"order\":[";for(long j=0;j<8;++j){if(j)std::cout<<',';std::cout<<d.order[j];}
  std::cout<<"],\"shapes\":[";for(long j=0;j<7;++j){if(j)std::cout<<',';std::cout<<d.shapes[j];}
  std::cout<<"],\"probability\":\""<<d.probability.str()<<"\",\"weight\":\""<<(d.weight*s.mass).str()<<"\"}";
 }
 std::cout<<']';
}
