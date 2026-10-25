// All order choices are summed exactly. No selector, probability divisor, or order average remains.
struct CompleteSharedOrder{std::array<long,6>order;long multiplicity;};
#include "shared_e3_full_sector_orbits_v1.hpp"
struct CompleteSharedMap{std::array<long,6>order;long multiplicity;std::array<long,5>shapes;std::array<Rat,6>alpha;Rat weight;};
std::vector<CompleteSharedMap> complete_shared_maps(const Graph&g,const std::array<Rat,5>&p,bool fourth){
 const auto&t=e3_mobius_table(g.id,fourth);std::vector<CompleteSharedMap>out;long total=0;
 for(const auto&orbit:shared_complete_orders(g.id)){
  CompleteSharedMap row;row.order=orbit.order;row.multiplicity=orbit.multiplicity;total+=row.multiplicity;
  long remaining=63;remaining^=1L<<row.order[0];row.alpha[row.order[0]]=Rat(1);Rat scale(1),jacobian(1);
  for(long j=1;j<6;++j){long shape=t.shapes[remaining];row.shapes[j-1]=shape;Rat ratio,derivative;
   if(shape<0){ratio=rat_power(p[j-1],-shape);derivative=Rat(-shape)*rat_power(p[j-1],-shape-1);}
   else{Rat denominator=Rat(1)+Rat(shape-1)*p[j-1];ratio=Rat(shape)*p[j-1]/denominator;derivative=Rat(shape)/rat_power(denominator,2);}
   scale=scale*ratio;row.alpha[row.order[j]]=scale;jacobian=jacobian*derivative*rat_power(ratio,5-j);remaining^=1L<<row.order[j];
  }
  if(remaining)throw std::runtime_error("complete six-edge order required");row.weight=Rat(row.multiplicity)*jacobian;out.push_back(row);
 }
 if(total!=720)throw std::runtime_error("all720 exact orders required");return out;
}
void complete_shared_map_output(const Graph&g,const std::array<Rat,5>&p,bool fourth){
 auto rows=complete_shared_maps(g,p,fourth);std::cout<<'[';bool first=true;
 for(const auto&row:rows){if(!first)std::cout<<',';first=false;std::cout<<"{\"order\":[";
  for(long j=0;j<6;++j){if(j)std::cout<<',';std::cout<<row.order[j];}
  std::cout<<"],\"multiplicity\":"<<row.multiplicity<<",\"shapes\":[";
  for(long j=0;j<5;++j){if(j)std::cout<<',';std::cout<<row.shapes[j];}
  std::cout<<"],\"alpha\":[";for(long j=0;j<6;++j){if(j)std::cout<<',';std::cout<<'\"'<<row.alpha[j].str()<<'\"';}
  std::cout<<"],\"weight\":\""<<row.weight.str()<<"\"}";
 }
 std::cout<<']';
}
std::map<Rat,Rat> prepare_complete_shared(const Graph&g,const std::array<Rat,5>&p,bool positive,bool fourth){
 if(!g.shared||g.states.at(0).weighted.size()!=1)throw std::invalid_argument("proved shared E3 graph required");
 std::map<Rat,Rat>arguments;const auto&weights=g.states.at(0).weighted;
 auto add=[&](const Rat&t,const Rat&c){if(!(Rat(0)<t))throw std::runtime_error("positive complete primitive argument required");
  auto it=arguments.find(t);if(it==arguments.end()){if(!fmpq_is_zero(c.v))arguments.emplace(t,c);}
  else{it->second=it->second+c;if(fmpq_is_zero(it->second.v))arguments.erase(it);}};
 for(const auto&row:complete_shared_maps(g,p,fourth)){
  Rat mass;for(const auto&a:row.alpha)mass=mass+a;Rat d=row.alpha[weights[0]],t0=(mass-d)/d,constant;
  for(const auto&state:g.states){if(state.weighted!=weights)throw std::runtime_error("shared edge identity required");
   auto uf=e3_exact_uf(state,row.alpha,positive);if(!(Rat(0)<uf[0])||!(Rat(0)<uf[1]))throw std::runtime_error("positive exact U,F required");
   Rat coefficient=row.weight*Rat(long(state.sign))/rat_power(uf[0],3);constant=constant+coefficient;
   add(t0+uf[1]/(uf[0]*d),Rat(0)-coefficient);
  }
  add(t0,constant);
 }
 Rat sum;for(const auto&entry:arguments)sum=sum+entry.second;
 if(!fmpq_is_zero(sum.v))throw std::runtime_error("primitive normalization must cancel exactly");return arguments;
}
#include "shared_square_small_native_v1.hpp"
#include "shared_square_cyclotomic_native_v1.hpp"
Real complete_shared_normalization(){
 static long cached_bits=0;static Real value;
 if(cached_bits!=precision){value=shared_square_P(Real(1),true);cached_bits=precision;}
 return value;
}
Real complete_shared_value(const Graph&g,const std::array<Rat,5>&p,bool positive,bool fourth,long requested_bits){
 const auto arguments=prepare_complete_shared(g,p,positive,fourth);std::string last;
 Rat small_coefficient;long small_count=0;
 for(const auto&entry:arguments)if(!(Rat("1/16")<entry.first)){small_coefficient=small_coefficient+entry.second;++small_count;}
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{Real total;long cyclic_count=0,closed_count=0;
   for(const auto&entry:arguments){Real t=to_real(entry.first),value;
    if(!(Rat("1/16")<entry.first))value=shared_square_small(t).value;
    else{Real w=4/sq(rt(t+4)+rt(t));auto evaluated=shared_square_cyclic(w);value=evaluated.value;
     if(evaluated.series)++cyclic_count;else ++closed_count;
    }
    total=total+to_real(entry.second)*value;
   }
   // A2(t)=P(w(t))-P(1); collect the normalization coefficient exactly once.
   if(!fmpq_is_zero(small_coefficient.v))total=total+to_real(small_coefficient)*complete_shared_normalization();
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){
    std::cerr<<"complete-shared-hybrid-working-bits "<<working<<" unique-primitives "<<arguments.size()
     <<" small-series "<<small_count<<" cyclic-series "<<cyclic_count<<" closed "<<closed_count
     <<" exact-normalization-coefficient "<<small_coefficient.str()<<'\n';return total;
   }
   last="complete hybrid shared forest missed80relativebits";
  }catch(const std::runtime_error&e){last=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("complete hybrid shared precision exhausted: "+last);
}
