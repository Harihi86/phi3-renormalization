// Four unresolved E2 graphs: full paired forest, analytic common scale.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <algorithm>
namespace shared_weighted_e2 {
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::Z;
using fourier_e2_root::readq;using fourier_e2_root::interval;
Real sq(const Real&a){Real r;arb_sqr(r.v,a.v,precision);return r;}
Real rt(const Real&a){Real r;arb_sqrt(r.v,a.v,precision);return r;}
Real lp(const Real&a){Real r;arb_log1p(r.v,a.v,precision);return r;}
Real real_log(const Real&a){Real r;arb_log(r.v,a.v,precision);return r;}
Real real_exp(const Real&a){Real r;arb_exp(r.v,a.v,precision);return r;}
Real em1(const Real&a){Real r;arb_expm1(r.v,a.v,precision);return r;}
Real power3(const Real&a){return a*sq(a);}
void pos(const Real&a,const char*msg){if(!arb_is_positive(a.v))throw std::runtime_error(msg);}
Real bubble_B(const Real&t){
 if(arb_is_zero(t.v))return Real(0);pos(t,"positive bubble argument");Real h;arb_asinh(h.v,(rt(t)/2).v,precision);
 return (t+4)*rt(t+4)/rt(t)*h/3-4*t/9-Real(4)/3;
}
struct State{int sign;std::array<unsigned,1> weighted;std::array<std::array<int,3>,5> routes;std::vector<std::array<int,2>> U,F;};
struct Graph{std::string id,sha;std::vector<State> states;};
#include "shared_weighted_e2_data.hpp"
Real det(const std::array<std::array<Real,3>,3>&a,unsigned n){
 std::array<unsigned,3> permutation={0,1,2};Real result;
 do{unsigned inversions=0;Real term(1);for(unsigned i=0;i<n;++i){term=term*a[i][permutation[i]];for(unsigned j=i+1;j<n;++j)inversions+=permutation[i]>permutation[j];}result=result+(inversions%2?-term:term);}
 while(std::next_permutation(permutation.begin(),permutation.begin()+n));return result;}
Real poly(const std::vector<std::array<int,2>>&terms,const std::array<Real,5>&alpha){
 Real total;for(const auto&term:terms){Real v(term[1]);for(unsigned j=0;j<5;++j)if(term[0]&(1<<j))v=v*alpha[j];total=total+v;}return total;}
std::array<Real,2> uf(const State&s,const std::array<Real,5>&alpha,bool positive){
 if(positive)return {poly(s.U,alpha),poly(s.F,alpha)};std::array<std::array<Real,3>,3>a;
 for(unsigned e=0;e<5;++e)for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)a[i][j]=a[i][j]+alpha[e]*Real(s.routes[e][i]*s.routes[e][j]);
 return {det(a,2),det(a,3)};}
Real density(const Graph&g,const std::array<Real,5>&alpha,const Real&v,bool positive){
 Real mass,total;for(const auto&a:alpha)mass=mass+a;auto weights=g.states.at(0).weighted;Real d=alpha[weights[0]];
 for(const auto&s:g.states){if(s.weighted!=weights)throw std::runtime_error("fixed shared edge required");auto pol=uf(s,alpha,positive);pos(pol[0],"positive U");
   Real z=pol[1]/pol[0],a=mass+v*z,b=bubble_B((a-d)/d);
   total=total+Real(s.sign)*sq(d)*sq(z)*(1-v)*sq(b)/(power3(a)*power3(pol[0]));}
 return total;}
Real value(const Graph&g,const std::array<Real,5>&p,bool positive,bool square){
 std::array<Real,4>x;Real jac(1),total;for(unsigned i=0;i<4;++i){x[i]=square?sq(p[i]):p[i];if(square)jac=jac*2*p[i];}
 for(unsigned pivot=0;pivot<5;++pivot){std::array<Real,5>alpha;unsigned next=0;for(unsigned e=0;e<5;++e)alpha[e]=e==pivot?Real(1):x[next++];total=total+density(g,alpha,p[4],positive);}
 total=total*jac;if(!arb_is_finite(total.v))throw std::runtime_error("finite full shared-bubble forest required");return total;}
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one double weighted E2 recipe required");
 auto g=load_weighted_e2(id);if(g.sha!=sha)throw std::invalid_argument("original double weighted identity required");return g;}
}
int main(int argc,char**argv){try{using namespace shared_weighted_e2;if(argc!=2)throw std::invalid_argument("one recipe required");auto g=recipe(argv[1]);
 std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||(mode!="shared-weighted-linear-original"&&mode!="shared-weighted-linear-positive"&&mode!="shared-weighted-square-original"&&mode!="shared-weighted-square-positive")||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded five-dimensional request required");
 precision=bits;bool positive=mode.find("positive")!=std::string::npos,square=mode.find("square")!=std::string::npos;
 std::cout<<"{\"id\":\""<<g.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(unsigned i=0;i<count;++i){std::array<Real,5>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open five-cube required");}if(i)std::cout<<',';interval(value(g,p,positive,square));}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing five-dimensional input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
