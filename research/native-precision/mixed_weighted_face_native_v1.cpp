// Complete weighted-bubble remainder; merge exact forest coefficients first.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <algorithm>
#include <map>
#include "fmpq.h"
namespace mixed_weighted_face {
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::interval;
#include "exact_rational_sector_native_v1.hpp"
#include "mixed_weighted_face_scalar_v1.hpp"
struct State{int sign;unsigned weighted;std::array<std::array<int,4>,8>routes;std::vector<std::array<int,2>>U,F;};
struct Graph{std::string id,sha;std::vector<State>states;};
#include "mixed_weighted_face_data.hpp"
#include "mixed_weighted_face_tables.hpp"
Rat polynomial(const std::vector<std::array<int,2>>&terms,const std::array<Rat,8>&a){
 Rat total;for(const auto&term:terms){Rat v{long(term[1])};for(long e=0;e<8;++e)if(term[0]&(1L<<e))v=v*a[e];total=total+v;}return total;
}
std::array<Rat,2> exact_uf(const State&s,const std::array<Rat,8>&a,bool positive){
 if(positive)return {polynomial(s.U,a),polynomial(s.F,a)};
 std::array<std::array<Rat,4>,4>m;
 for(long e=0;e<8;++e)for(long i=0;i<4;++i)for(long j=0;j<4;++j)m[i][j]=m[i][j]+a[e]*Rat(long(s.routes[e][i]*s.routes[e][j]));
 Rat determinant(1);
 for(long k=0;k<3;++k){Rat pivot=m[k][k];if(!(Rat(0)<pivot))throw std::runtime_error("positive Gaussian pivot required");determinant=determinant*pivot;
  for(long i=k+1;i<4;++i)for(long j=k+1;j<4;++j)m[i][j]=m[i][j]-m[i][k]*m[k][j]/pivot;
 }
 return {determinant,determinant*m[3][3]};
}
struct Term{long kind;Rat argument,coefficient;};
std::vector<Term> prepare(const Graph&g,const std::array<Rat,8>&point,bool positive,bool third){
 const auto&table=projective_table(g.id,third);std::vector<Term>terms;
 for(long orbit=0;orbit<8;++orbit){Rat selector=point[0]+Rat(orbit)/Rat(8);if(!(selector<Rat(1)))selector=selector-Rat(1);
  auto mapped=rational_sector_point(table,selector,point);Rat mass;for(const auto&a:mapped.alpha)mass=mass+a;
  unsigned weighted=g.states.at(0).weighted;Rat d=mapped.alpha[weighted],t0=(mass-d)/d;
  std::map<std::pair<long,Rat>,Rat>combined;
  auto add=[&](long kind,const Rat&t,const Rat&co){auto key=std::make_pair(kind,t);combined[key]=combined[key]+co;};
  for(const auto&state:g.states){
   if(state.weighted!=weighted)throw std::runtime_error("original weighted edge identity required");
   auto v=exact_uf(state,mapped.alpha,positive);
   if(!(Rat(0)<v[0])||!(Rat(0)<v[1]))throw std::runtime_error("positive exact U and F required");
   Rat t1=t0+v[1]/(d*v[0]),co=Rat(long(state.sign))*d/rat_power(v[0],3);
   add(0,t1,co*t1);add(1,t1,Rat(-1)*co);add(0,t0,Rat(-1)*co*t1);add(1,t0,co);
  }
  for(const auto&row:combined)if(!fmpq_is_zero(row.second.v))terms.push_back({row.first.first,row.first.second,row.second*mapped.weight/Rat(8)});
 }
 return terms;
}
Real value(const std::vector<Term>&terms,long bits){
 for(long work=bits;;work=std::min(work*2,8192L)){
  precision=work;Real total;
  for(const auto&term:terms){Real t=to_real(term.argument);total=total+to_real(term.coefficient)*(term.kind==0?bubble_A(t):weighted_C(t));}
  if(arb_is_finite(total.v)&&(arb_is_zero(total.v)||arb_rel_accuracy_bits(total.v)>=80)){std::cerr<<"mixed-face-working-bits "<<work<<'\n';return total;}
  if(work==8192)break;
 }
 throw std::runtime_error("complete weighted forest failed bounded 80-relative-bit arithmetic gate");
}
void output_coefficients(const std::vector<Term>&terms){std::cout<<'[';bool first=true;
 for(const auto&t:terms){if(!first)std::cout<<',';first=false;std::cout<<'['<<t.kind<<",\""<<t.argument.str()<<"\",\""<<t.coefficient.str()<<"\"]";}std::cout<<']';}
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original mixed-graph identity required");
 auto g=load_projective(id);if(g.sha!=sha)throw std::invalid_argument("original graph hash required");return g;}
}
int main(int argc,char**argv){try{using namespace mixed_weighted_face;
 if(argc!=2)throw std::invalid_argument("one recipe required");auto graph=recipe(argv[1]);std::string mode;long bits,count;
 if(!(std::cin>>mode>>bits>>count)||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded exact request required");
 const std::vector<std::string>modes={"face2-original","face2-positive","face3-original","face3-positive","map-face2","map-face3","coeff-face2-original","coeff-face2-positive","coeff-face3-original","coeff-face3-positive"};
 if(std::find(modes.begin(),modes.end(),mode)==modes.end())throw std::invalid_argument("declared mixed-face mode required");
 bool third=mode.find("face3")!=std::string::npos,positive=mode.find("positive")!=std::string::npos,maps=mode.find("map-")==0,coeff=mode.find("coeff-")==0;
 std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\""<<(maps?"maps":coeff?"coefficients":"values")<<"\":[";
 for(long i=0;i<count;++i){std::array<Rat,8>point;for(auto&v:point){v=read_rational();if(!(Rat(0)<v)||!(v<Rat(1)))throw std::invalid_argument("strict open exact 8-cube required");}
  if(i)std::cout<<',';precision=bits;
  if(maps)rational_sector_map_output(projective_table(graph.id,third),point);
  else{auto terms=prepare(graph,point,positive,third);if(coeff)output_coefficients(terms);else interval(value(terms,bits));}}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
