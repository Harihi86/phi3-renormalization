// Exact conditional partitions and Jacobians; included inside a numerical namespace.
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
template<std::size_t N>struct RationalSectorTable{std::array<long,(1UL<<N)>shapes,weights;long normalization,scale;};
template<std::size_t N>struct RationalSectorPoint{std::array<long,N>order;std::array<long,N-1>shapes;std::array<Rat,N>alpha;Rat probability{1},weight;};
template<std::size_t N>RationalSectorPoint<N> rational_sector_point(const RationalSectorTable<N>&t,Rat selector,const std::array<Rat,N>&p){
 if(selector<Rat(0)||!(selector<Rat(1)))throw std::invalid_argument("half-open exact selector required");
 RationalSectorPoint<N>out;long remaining=(1L<<N)-1;
 for(long depth=0;depth<long(N);++depth){
  long total=0;for(long edge=0;edge<long(N);++edge)if(remaining&(1L<<edge))total+=t.weights[remaining^(1L<<edge)];
  Rat left;bool found=false;
  for(long edge=0;edge<long(N);++edge)if(remaining&(1L<<edge)){
   Rat probability=Rat(t.weights[remaining^(1L<<edge)])/Rat(total);
   if(selector<left+probability){selector=(selector-left)/probability;out.probability=out.probability*probability;out.order[depth]=edge;remaining^=1L<<edge;if(remaining)out.shapes[depth]=t.shapes[remaining];found=true;break;}
   left=left+probability;
  }
  if(!found)throw std::runtime_error("exact conditional CDF gap");
 }
 out.alpha[out.order[0]]=Rat(1);Rat scale(1),jacobian(1);
 for(long j=1;j<long(N);++j){
  long shape=out.shapes[j-1];Rat ratio,derivative;
  if(shape<0){long power=-shape;ratio=rat_power(p[j],power);derivative=Rat(power)*rat_power(p[j],power-1);}
  else{Rat denominator=Rat(1)+Rat(shape-1)*p[j];ratio=Rat(shape)*p[j]/denominator;derivative=Rat(shape)/rat_power(denominator,2);}
  scale=scale*ratio;out.alpha[out.order[j]]=scale;jacobian=jacobian*derivative*rat_power(ratio,long(N)-j-1);
 }
 out.weight=jacobian/out.probability;return out;
}
template<std::size_t N>Rat rational_determinant(const std::array<std::array<Rat,N>,N>&m){
 std::array<long,N>perm;for(long i=0;i<long(N);++i)perm[i]=i;Rat total;
 do{Rat term(1);long inversions=0;for(long i=0;i<long(N);++i){term=term*m[i][perm[i]];for(long j=i+1;j<long(N);++j)inversions+=perm[i]>perm[j];}total=total+Rat(inversions%2?-1L:1L)*term;}
 while(std::next_permutation(perm.begin(),perm.end()));return total;
}
template<std::size_t N>void rational_sector_map_output(const RationalSectorTable<N>&t,const std::array<Rat,N>&p){
 std::cout<<'[';
 for(long orbit=0;orbit<long(N);++orbit){Rat selector=p[0]+Rat(orbit)/Rat(long(N));if(!(selector<Rat(1)))selector=selector-Rat(1);auto d=rational_sector_point(t,selector,p);
  if(orbit)std::cout<<',';std::cout<<"{\"alpha\":[";for(long j=0;j<long(N);++j){if(j)std::cout<<',';std::cout<<'\"'<<d.alpha[j].str()<<'\"';}
  std::cout<<"],\"order\":[";for(long j=0;j<long(N);++j){if(j)std::cout<<',';std::cout<<d.order[j];}
  std::cout<<"],\"shapes\":[";for(long j=0;j<long(N)-1;++j){if(j)std::cout<<',';std::cout<<d.shapes[j];}
  std::cout<<"],\"probability\":\""<<d.probability.str()<<"\",\"weight\":\""<<d.weight.str()<<"\"}";
 }
 std::cout<<']';
}
