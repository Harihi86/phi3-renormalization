// Exact rational discrete sectors, exact independent U/F evaluations, bounded precision recovery.
struct Rat{
 fmpq_t v;Rat(){fmpq_init(v);}explicit Rat(long x){fmpq_init(v);fmpq_set_si(v,x,1);}
 explicit Rat(const std::string&s){fmpq_init(v);if(fmpq_set_str(v,s.c_str(),10)||fmpz_sgn(fmpq_denref(v))<=0){fmpq_clear(v);throw std::invalid_argument("rational with positive denominator required");}fmpq_canonicalise(v);}
 Rat(const Rat&o){fmpq_init(v);fmpq_set(v,o.v);}Rat&operator=(const Rat&o){fmpq_set(v,o.v);return *this;}~Rat(){fmpq_clear(v);}
 std::string str()const{char*p=fmpq_get_str(nullptr,10,v);std::string s(p);flint_free(p);return s;}
};
Rat operator+(const Rat&a,const Rat&b){Rat r;fmpq_add(r.v,a.v,b.v);return r;}
Rat operator-(const Rat&a,const Rat&b){Rat r;fmpq_sub(r.v,a.v,b.v);return r;}
Rat operator*(const Rat&a,const Rat&b){Rat r;fmpq_mul(r.v,a.v,b.v);return r;}
Rat operator/(const Rat&a,const Rat&b){if(fmpq_is_zero(b.v))throw std::runtime_error("nonzero rational divisor required");Rat r;fmpq_div(r.v,a.v,b.v);return r;}
bool operator<(const Rat&a,const Rat&b){return fmpq_cmp(a.v,b.v)<0;}
Real to_real(const Rat&a){Real r;arb_set_fmpq(r.v,a.v,precision);return r;}
Rat rat_power(Rat a,long n){if(n<0)throw std::invalid_argument("nonnegative exact power required");Rat r(1);while(n){if(n&1)r=r*a;n>>=1;if(n)a=a*a;}return r;}
Rat read_rational(){std::string n,d;if(!(std::cin>>n>>d)||n.size()>10000||d.size()>10000)throw std::invalid_argument("bounded exact rational required");return Rat(n+"/"+d);}
struct HeppTable{std::array<long,256> powers,weights;long normalization;};
#include "mixed_forest_hepp_tables_v1.hpp"
struct HeppPoint{std::array<long,8>order;std::array<long,7>powers;std::array<Rat,8>alpha;Rat probability{1},weight;};
HeppPoint hepp_point(const HeppTable&t,Rat selector,const std::array<Rat,8>&p){
 if(selector<Rat(0)||!(selector<Rat(1)))throw std::invalid_argument("half-open exact selector required");
 HeppPoint out;long remaining=255;
 for(long depth=0;depth<8;++depth){
  long total=0;for(long edge=0;edge<8;++edge)if(remaining&(1L<<edge))total+=t.weights[remaining^(1L<<edge)];
  Rat left;bool found=false;
  for(long edge=0;edge<8;++edge)if(remaining&(1L<<edge)){
   Rat probability=Rat(t.weights[remaining^(1L<<edge)])/Rat(total);
   if(selector<left+probability){selector=(selector-left)/probability;out.probability=out.probability*probability;out.order[depth]=edge;remaining^=1L<<edge;if(remaining)out.powers[depth]=t.powers[remaining];found=true;break;}
   left=left+probability;
  }
  if(!found)throw std::runtime_error("exact Hepp conditional CDF gap");
 }
 out.alpha[out.order[0]]=Rat(1);Rat scale(1);out.weight=Rat(t.normalization);
 for(long j=1;j<8;++j){Rat ratio=rat_power(p[j],out.powers[j-1]);scale=scale*ratio;out.alpha[out.order[j]]=scale;out.weight=out.weight*rat_power(ratio,8-j)/p[j];}
 return out;
}
Rat rat_determinant(const std::array<std::array<Rat,4>,4>&m,long size){
 std::array<long,4>perm={0,1,2,3};Rat total;
 do{Rat term(1);long inversions=0;for(long i=0;i<size;++i){term=term*m[i][perm[i]];for(long j=i+1;j<size;++j)inversions+=perm[i]>perm[j];}total=total+Rat(inversions%2?-1L:1L)*term;}
 while(std::next_permutation(perm.begin(),perm.begin()+size));return total;
}
Rat rat_polynomial(const std::vector<std::array<int,2>>&terms,const std::array<Rat,8>&a){Rat sum;for(const auto&term:terms){Rat product{long(term[1])};for(long edge=0;edge<8;++edge)if(term[0]&(1L<<edge))product=product*a[edge];sum=sum+product;}return sum;}
std::array<Rat,2> exact_uf(const State&s,const std::array<Rat,8>&a,bool positive){
 if(positive)return {rat_polynomial(s.U,a),rat_polynomial(s.F,a)};
 std::array<std::array<Rat,4>,4>m;
 for(long edge=0;edge<8;++edge)for(long i=0;i<4;++i)for(long j=0;j<4;++j)m[i][j]=m[i][j]+a[edge]*Rat(long(s.routes[edge][i]*s.routes[edge][j]));
 return {rat_determinant(m,3),rat_determinant(m,4)};
}
struct HeppTerm{Rat t1,coefficient;};
struct HeppSector{Rat t0,weight;std::vector<HeppTerm>terms;};
std::vector<HeppSector> prepare_hepp(const Graph&g,const std::array<Rat,8>&p,bool positive,bool fourth){
 const auto&t=hepp_table(g.id,fourth);std::vector<HeppSector>sectors;
 for(long orbit=0;orbit<8;++orbit){Rat selector=p[0]+Rat(orbit)/Rat(8);if(!(selector<Rat(1)))selector=selector-Rat(1);auto decoded=hepp_point(t,selector,p);
  Rat mass;for(const auto&a:decoded.alpha)mass=mass+a;Rat d=decoded.alpha[g.states.at(0).weighted];HeppSector sector;sector.t0=(mass-d)/d;sector.weight=decoded.weight/Rat(8);
  for(const auto&state:g.states){auto uf=exact_uf(state,decoded.alpha,positive);if(!(Rat(0)<uf[0])||!(Rat(0)<uf[1]))throw std::runtime_error("positive exact U and F required");
   Rat delta=uf[1]/(uf[0]*d);sector.terms.push_back({sector.t0+delta,Rat(long(state.sign))*d/rat_power(uf[0],3)});}
  sectors.push_back(sector);
 }
 return sectors;
}
Real hepp_value(const Graph&g,const std::array<Rat,8>&p,bool positive,bool fourth,long requested_bits){
 auto sectors=prepare_hepp(g,p,positive,fourth);std::string last_error;
 for(long working=requested_bits;;working=std::min(working*2,8192L)){
  precision=working;
  try{
   Real total;
   for(const auto&sector:sectors){Real t0=to_real(sector.t0),a0=bubble_A(t0),c0=weighted_C(t0),subtotal;
    for(const auto&term:sector.terms){Real t1=to_real(term.t1);subtotal=subtotal+to_real(term.coefficient)*(t1*(bubble_A(t1)-a0)-(weighted_C(t1)-c0));}
    total=total+to_real(sector.weight)*subtotal;
   }
   if(arb_is_finite(total.v)&&arb_rel_accuracy_bits(total.v)>=80){std::cerr<<"hepp-working-bits "<<working<<'\n';return total;}
   last_error="complete Hepp remainder missed 80 relative bits";
  }catch(const std::runtime_error&e){last_error=e.what();}
  if(working==8192)break;
 }
 throw std::runtime_error("bounded Hepp precision recovery exhausted: "+last_error);
}
void hepp_map_output(const Graph&g,const std::array<Rat,8>&p,bool fourth){
 const auto&t=hepp_table(g.id,fourth);std::cout<<'[';
 for(long orbit=0;orbit<8;++orbit){Rat selector=p[0]+Rat(orbit)/Rat(8);if(!(selector<Rat(1)))selector=selector-Rat(1);auto d=hepp_point(t,selector,p);
  if(orbit)std::cout<<',';std::cout<<"{\"alpha\":[";for(long j=0;j<8;++j){if(j)std::cout<<',';std::cout<<'\"'<<d.alpha[j].str()<<'\"';}
  std::cout<<"],\"order\":[";for(long j=0;j<8;++j){if(j)std::cout<<',';std::cout<<d.order[j];}
  std::cout<<"],\"powers\":[";for(long j=0;j<7;++j){if(j)std::cout<<',';std::cout<<d.powers[j];}
  std::cout<<"],\"probability\":\""<<d.probability.str()<<"\",\"weight\":\""<<d.weight.str()<<"\"}";
 }
 std::cout<<']';
}
