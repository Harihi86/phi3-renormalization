// Complete common-scale remainder with a proved rho=2 ellipse enclosure.
struct CCRule{long degree;std::vector<Real> nodes,weights;Real error_factor;};
Real cc_abs(const Real&a){Real r;arb_abs(r.v,a.v);return r;}
Real cc_symmetric(const Real&a){Real r;arb_union(r.v,(-a).v,a.v,precision);return r;}
const CCRule& cc_rule(long n){
 static std::vector<CCRule> cache;
 for(const auto&r:cache)if(r.degree==n)return r;
 if(n!=64&&n!=128&&n!=256&&n!=512)throw std::invalid_argument("predeclared CC degree required");
 CCRule r;r.degree=n;std::vector<Real> cosines;
 for(long j=0;j<2*n;++j){Real c,a=Real(j)/Real(n);arb_cos_pi(c.v,a.v,precision);cosines.push_back(c);}
 for(long j=0;j<=n;++j){
  Real w(1);for(long k=1;k<n/2;++k)w=w-2*cosines[(2*k*j)%(2*n)]/Real(4*k*k-1);
  w=w-cosines[(n*j)%(2*n)]/Real(n*n-1);w=w*Real(j==0||j==n?1L:2L)/Real(2*n);
  pos(w,"positive normalized CC weight");r.nodes.push_back((1+cosines[j])/2);r.weights.push_back(w);
 }
 r.error_factor=Real(3468)/7595;arb_mul_2exp_si(r.error_factor.v,r.error_factor.v,-n);
 cache.push_back(r);return cache.back();
}
Real cc_weight(const Real&t){
 pos(t,"strict positive scale argument");Real h;arb_asinh(h.v,(rt(t)/2).v,precision);
 return ((t+4)*rt((t+4)/t)*h/3-4*t/9-Real(4)/3)/sq(1+t);
}
struct CCScaleTerm{Real t0,delta,coefficient;};
Real cc_value(const Graph&g,bool positive,bool square,const std::array<Real,7>&p){
 std::array<Real,7>x;Real jacobian(1);std::vector<CCScaleTerm> terms;Real error_scale;
 for(long j=0;j<7;++j){x[j]=square?sq(p[j]):p[j];if(square)jacobian=jacobian*2*p[j];}
 for(long pivot=0;pivot<8;++pivot){
  std::array<Real,8>alpha;long next=0;Real mass;
  for(long e=0;e<8;++e){alpha[e]=e==pivot?Real(1):x[next++];mass=mass+alpha[e];}
  unsigned weighted=g.states.at(0).weighted;Real edge=alpha[weighted],t0=(mass-edge)/edge;
  for(const auto&state:g.states){
   if(state.weighted!=weighted)throw std::runtime_error("same weighted edge required");
   auto pair=uf(state,alpha,positive);pos(pair[0],"positive U");pos(pair[1],"positive F");
   Real q=pair[1]/pair[0],delta=q/edge,coefficient=Real(state.sign)*edge*sq(delta)/power3(pair[0]);
   terms.push_back({t0,delta,coefficient});error_scale=error_scale+cc_abs(coefficient);
  }
 }
 for(long n:{64L,128L,256L,512L}){
  const auto&r=cc_rule(n);Real total;
  for(long j=0;j<=n;++j){Real sum;for(const auto&t:terms)sum=sum+t.coefficient*cc_weight(t.t0+t.delta*r.nodes[j]);
   total=total+r.weights[j]*(1-r.nodes[j])*sum;}
  Real value=(total+cc_symmetric(error_scale*r.error_factor))*jacobian;
  if(!arb_is_finite(value.v))throw std::runtime_error("finite certified common-scale remainder required");
  if(arb_rel_accuracy_bits(value.v)>=80){std::cerr<<"cc-scale-degree "<<n<<'\n';return value;}
 }
 throw std::runtime_error("predeclared complete-scale enclosure missed 80 relative bits");
}
void cc_rule_moments(){
 std::cout<<"{\"rules\":[";bool first=true;
 for(long n:{64L,128L,256L,512L}){
  const auto&r=cc_rule(n);std::vector<Real> moments(n+1);
  for(long j=0;j<=n;++j){Real x=2*r.nodes[j]-1,t0(1),t1=x;moments[0]=moments[0]+r.weights[j];moments[1]=moments[1]+r.weights[j]*t1;
   for(long k=2;k<=n;++k){Real t2=2*x*t1-t0;moments[k]=moments[k]+r.weights[j]*t2;t0=t1;t1=t2;}}
  if(!first)std::cout<<',';first=false;std::cout<<"{\"degree\":"<<n<<",\"moments\":[";
  for(long k=0;k<=n;++k){if(k)std::cout<<',';interval(moments[k]);}std::cout<<"],\"error_factor\":";interval(r.error_factor);std::cout<<'}';
 }
 std::cout<<"]}\n";
}
