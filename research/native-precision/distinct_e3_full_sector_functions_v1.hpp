// All order choices are summed exactly. No selector, probability divisor, or order average remains.
struct CompleteDistinctOrder{std::array<long,6>order;long multiplicity;};
#include "distinct_e3_full_sector_orbits_v1.hpp"
struct CompleteDistinctMap{std::array<long,6>order;long multiplicity;std::array<long,5>shapes;std::array<Rat,6>alpha;Rat weight;};
std::vector<CompleteDistinctMap> complete_distinct_maps(const Graph&g,const std::array<Rat,5>&p,bool fourth){
 const auto&t=e3_mobius_table(g.id,fourth);std::vector<CompleteDistinctMap>out;long total=0;
 for(const auto&orbit:distinct_complete_orders(g.id)){
  CompleteDistinctMap row;row.order=orbit.order;row.multiplicity=orbit.multiplicity;total+=row.multiplicity;
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
void complete_distinct_map_output(const Graph&g,const std::array<Rat,5>&p,bool fourth){
 auto rows=complete_distinct_maps(g,p,fourth);std::cout<<'[';bool first=true;
 for(const auto&row:rows){if(!first)std::cout<<',';first=false;std::cout<<"{\"order\":[";
  for(long j=0;j<6;++j){if(j)std::cout<<',';std::cout<<row.order[j];}
  std::cout<<"],\"multiplicity\":"<<row.multiplicity<<",\"shapes\":[";
  for(long j=0;j<5;++j){if(j)std::cout<<',';std::cout<<row.shapes[j];}
  std::cout<<"],\"alpha\":[";for(long j=0;j<6;++j){if(j)std::cout<<',';std::cout<<'\"'<<row.alpha[j].str()<<'\"';}
  std::cout<<"],\"weight\":\""<<row.weight.str()<<"\"}";
 }
 std::cout<<']';
}

using DistinctArguments=std::map<Rat,std::array<Rat,2>>;
DistinctArguments prepare_complete_distinct(const Graph&g,const std::array<Rat,6>&p,bool positive,bool fourth){
 if(g.shared||g.states.at(0).weighted.size()!=2)throw std::invalid_argument("two ordered distinct weighted edges required");
 const auto&weights=g.states.at(0).weighted;if(weights[0]==weights[1])throw std::invalid_argument("distinct edge roles required");
 std::array<Rat,5>ratios;for(long j=0;j<5;++j)ratios[j]=p[j];
 Rat v=p[5],r=Rat(1)-v*v,tau=(Rat(1)-r*r)/Rat(4),spectral=v*rat_power(r,4)/Rat(6);
 Rat c=(tau-Rat(2))/rat_power(Rat(1)-tau,2),h=Rat(1)/(Rat(1)-tau),k=Rat(1)/(tau*rat_power(Rat(1)-tau,2));
 DistinctArguments arguments;
 auto add=[&](const Rat&t,long kind,const Rat&coefficient){
  if(!(Rat(0)<t))throw std::runtime_error("positive distinct primitive argument required");
  auto&entry=arguments[t];entry[kind]=entry[kind]+coefficient;
  if(fmpq_is_zero(entry[0].v)&&fmpq_is_zero(entry[1].v))arguments.erase(t);
 };
 for(const auto&row:complete_distinct_maps(g,ratios,fourth)){
  Rat mass;for(const auto&a:row.alpha)mass=mass+a;
  Rat d=row.alpha[weights[0]],e=row.alpha[weights[1]],t0=(mass-d)/d,s0=(mass+e*(Rat(1)/tau-Rat(1))-d)/d;
  for(const auto&state:g.states){
   if(state.weighted!=weights)throw std::runtime_error("ordered weighted-edge identities changed");
   auto uf=e3_exact_uf(state,row.alpha,positive);
   if(!(Rat(0)<uf[0])||!(Rat(0)<uf[1]))throw std::runtime_error("positive exact distinct U,F required");
   Rat shift=uf[1]/(uf[0]*d),t1=t0+shift,s1=s0+shift;
   Rat scale=row.weight*spectral*Rat(long(state.sign))/rat_power(uf[0],3);
   add(t0,0,scale*c);add(t1,0,Rat(0)-scale*c);add(s0,0,scale*k);add(s1,0,Rat(0)-scale*k);
   Rat derivative=Rat(0)-scale*h*e/d;
   add(t0,1,derivative/rat_power(Rat(1)+t0,2));add(t1,1,(Rat(0)-derivative)/rat_power(Rat(1)+t1,2));
  }
 }
 Rat sum;for(const auto&entry:arguments)sum=sum+entry.second[0];
 if(!fmpq_is_zero(sum.v))throw std::runtime_error("distinct primitive normalization must cancel exactly");
 return arguments;
}
void complete_distinct_coeff_output(const Graph&g,const std::array<Rat,6>&p,bool positive,bool fourth){
 const auto arguments=prepare_complete_distinct(g,p,positive,fourth);std::cout<<'[';bool first=true;
 for(const auto&entry:arguments){if(!first)std::cout<<',';first=false;
  std::cout<<"{\"t\":\""<<entry.first.str()<<"\",\"A\":\""<<entry.second[0].str()<<"\",\"B\":\""<<entry.second[1].str()<<"\"}";}
 std::cout<<']';
}
Real complete_distinct_value(const Graph&g,const std::array<Rat,6>&p,bool positive,bool fourth,long requested_bits){
 const auto arguments=prepare_complete_distinct(g,p,positive,fourth);std::string last;
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{Real total;
   for(const auto&entry:arguments){Real t=to_real(entry.first);
    if(!fmpq_is_zero(entry.second[0].v))total=total+to_real(entry.second[0])*bubble_A(t);
    if(!fmpq_is_zero(entry.second[1].v))total=total+to_real(entry.second[1])*bubble_B(t);
   }
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){std::cerr<<"complete-distinct-working-bits "<<working<<" unique-arguments "<<arguments.size()<<'\n';return total;}
   last="complete distinct forest missed80relativebits";
  }catch(const std::runtime_error&e){last=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("complete distinct precision exhausted: "+last);
}
