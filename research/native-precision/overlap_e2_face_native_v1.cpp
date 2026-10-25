// Original overlapping T0 forests with the complete two-point T2 remainder.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <algorithm>
#include "fmpq.h"
namespace overlap_e2_face {
using arb_cw::Real;using arb_cw::precision;using fourier_e2_root::interval;
#include "exact_rational_sector_native_v1.hpp"
struct State{int sign;std::array<std::array<int,5>,11>routes;std::vector<std::array<int,2>>U,F;};
struct Graph{std::string id,sha;std::vector<State>states;std::vector<std::array<int,2>>pairs;std::vector<int>unpaired;};
#include "overlap_e2_face_data.hpp"
#include "overlap_e2_face_tables.hpp"
Rat polynomial(const std::vector<std::array<int,2>>&terms,const std::array<Rat,11>&a){
 Rat total;for(const auto&term:terms){Rat value{long(term[1])};for(long e=0;e<11;++e)if(term[0]&(1L<<e))value=value*a[e];total=total+value;}return total;
}
std::array<Rat,2> exact_uf(const State&s,const std::array<Rat,11>&a,bool positive){
 if(positive)return {polynomial(s.U,a),polynomial(s.F,a)};
 std::array<std::array<Rat,5>,5>m;
 for(long e=0;e<11;++e)for(long i=0;i<5;++i)for(long j=0;j<5;++j)m[i][j]=m[i][j]+a[e]*Rat(long(s.routes[e][i]*s.routes[e][j]));
 Rat determinant(1);
 for(long k=0;k<4;++k){Rat pivot=m[k][k];if(!(Rat(0)<pivot))throw std::runtime_error("positive Gaussian pivot required");determinant=determinant*pivot;
  for(long i=k+1;i<5;++i)for(long j=k+1;j<5;++j)m[i][j]=m[i][j]-m[i][k]*m[k][j]/pivot;
 }
 return {determinant,determinant*m[4][4]};
}
struct Term{int kind;Rat argument,coefficient;}; // kind 0: phi, kind 1: log1p
std::vector<Term> prepare(const Graph&g,const std::array<Rat,11>&point,bool positive,bool fourth){
 const auto&t=projective_table(g.id,fourth);std::vector<Term>terms;
 for(long orbit=0;orbit<11;++orbit){Rat selector=point[0]+Rat(orbit)/Rat(11);if(!(selector<Rat(1)))selector=selector-Rat(1);
  auto mapped=rational_sector_point(t,selector,point);Rat mass;for(const auto&a:mapped.alpha)mass=mass+a;
  std::vector<std::array<Rat,2>>values;
  for(const auto&state:g.states){auto v=exact_uf(state,mapped.alpha,positive);
   if(!(Rat(0)<v[0])||!(Rat(0)<v[1]))throw std::runtime_error("positive exact U and F required");values.push_back({v[0],v[1]/(mass*v[0])});}
  std::vector<int>seen(g.states.size(),0);
  for(const auto&pair:g.pairs){const auto&b=values.at(pair[0]);const auto&z=values.at(pair[1]);
   if(g.states.at(pair[0]).sign!=-g.states.at(pair[1]).sign)throw std::runtime_error("opposite paired forest signs required");
   ++seen.at(pair[0]);++seen.at(pair[1]);
   Rat scale=mapped.weight*mass*Rat(long(g.states.at(pair[0]).sign))/Rat(11),ib=Rat(1)/rat_power(b[0],3),iz=Rat(1)/rat_power(z[0],3);
   Rat delta=b[1]-z[1],relative=delta/(Rat(1)+z[1]);if(!(Rat(-1)<relative))throw std::runtime_error("positive relative log domain required");
   terms.push_back({0,z[1],scale*(ib-iz)});
   terms.push_back({1,z[1],scale*ib*delta});
   terms.push_back({0,relative,scale*ib*(Rat(1)+z[1])});
  }
  for(int index:g.unpaired){const auto&v=values.at(index);++seen.at(index);
   terms.push_back({0,v[1],mapped.weight*mass*Rat(long(g.states.at(index).sign))/(Rat(11)*rat_power(v[0],3))});}
  if(std::any_of(seen.begin(),seen.end(),[](int n){return n!=1;}))throw std::runtime_error("every original forest required exactly once");
 }
 return terms;
}
Real value(const std::vector<Term>&terms,long bits){
 for(long work=bits;;work=std::min(work*2,8192L)){
  precision=work;Real total;
  for(const auto&term:terms){if(fmpq_is_zero(term.coefficient.v)||fmpq_is_zero(term.argument.v))continue;
   Real argument=to_real(term.argument),function;arb_log1p(function.v,argument.v,precision);
   if(term.kind==0)function=(Real(1L)+argument)*function-argument;
   else if(term.kind!=1)throw std::runtime_error("unknown scalar remainder");
   total=total+to_real(term.coefficient)*function;
  }
  if(arb_is_finite(total.v)&&(arb_is_zero(total.v)||arb_rel_accuracy_bits(total.v)>=80)){std::cerr<<"projective-working-bits "<<work<<'\n';return total;}
  if(work==8192)break;
 }
 throw std::runtime_error("complete paired forest failed bounded 80-relative-bit arithmetic gate");
}
void output_coefficients(const std::vector<Term>&terms){std::cout<<'[';bool first=true;
 for(const auto&t:terms){if(!first)std::cout<<',';first=false;std::cout<<"["<<t.kind<<",\""<<t.argument.str()<<"\",\""<<t.coefficient.str()<<"\"]";}std::cout<<']';}
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;
 if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original two-point identity required");
 auto g=load_projective(id);if(g.sha!=sha)throw std::invalid_argument("original graph hash required");return g;}
}
int main(int argc,char**argv){try{using namespace overlap_e2_face;
 if(argc!=2)throw std::invalid_argument("one recipe required");auto graph=recipe(argv[1]);std::string mode;long bits,count;
 if(!(std::cin>>mode>>bits>>count)||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded exact request required");
 const std::vector<std::string>modes={"face2-original","face2-positive","face3-original","face3-positive","map-face2","map-face3","coeff-face2-original","coeff-face2-positive","coeff-face3-original","coeff-face3-positive"};
 if(std::find(modes.begin(),modes.end(),mode)==modes.end())throw std::invalid_argument("declared projective mode required");
 bool fourth=mode.find("face3")!=std::string::npos,positive=mode.find("positive")!=std::string::npos,maps=mode.find("map-")==0,coeff=mode.find("coeff-")==0;
 std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\""<<(maps?"maps":coeff?"coefficients":"values")<<"\":[";
 for(long i=0;i<count;++i){std::array<Rat,11>point;for(auto&v:point){v=read_rational();if(!(Rat(0)<v)||!(v<Rat(1)))throw std::invalid_argument("strict open exact 11-cube required");}
  if(i)std::cout<<',';precision=bits;
  if(maps)rational_sector_map_output(projective_table(graph.id,fourth),point);
  else{auto terms=prepare(graph,point,positive,fourth);if(coeff)output_coefficients(terms);else interval(value(terms,bits));}}
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
