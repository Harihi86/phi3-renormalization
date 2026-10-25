// Six mixed unresolved E2 graphs: full paired forest, analytic common scale.
#include "disjoint-control-residual-native-v2/frozen_moments.hpp"
#include <algorithm>
namespace weighted_e2 {
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
Z clog(const Z&a){Z r;acb_log(r.v,a.v,precision);return r;}
Z dilog(const Z&a){Z r;acb_polylog_si(r.v,2,a.v,precision);return r;}
Real atan_log(const Real&q){
    Z a;arb_sqrt_ui(acb_imagref(a.v),3,precision);Z sum;
    for(long b:{-1L,1L}){Z zz=-a/(Z(b)-a);Real u=1+Real(b)*q;pos(u,"positive marginal logarithm argument");
        sum=sum+(clog(Z(1)-a/Z(b))*Z(real_log(u))-dilog(zz*Z(u))+dilog(zz))/Z(b);}
    Real r;arb_set(r.v,acb_imagref(sum.v));return r/2;
}
Real bubble_A(const Real&t){
    if(arb_is_zero(t.v))return Real(0);pos(t,"nonnegative bubble primitive argument");
    Real q=rt(t/(t+4)),h,angle;pos(1-q,"bubble primitive finite q");
    arb_asinh(h.v,(rt(t)/2).v,precision);arb_atan(angle.v,(rt(Real(3))*q).v,precision);
    Real value=sq(h)/3+2*rt(Real(3))*(h*angle-atan_log(q))+4*q*h/(1+3*sq(q))-Real(17)/18*lp(t)-Real(32)/9*sq(q)/(1+3*sq(q));
    if(!arb_is_finite(value.v))throw std::runtime_error("finite bubble primitive required");return value;
}

Real weighted_C(const Real&t){
    if(arb_is_zero(t.v))return Real(0);pos(t,"nonnegative weighted primitive argument");
    Real h;arb_asinh(h.v,(rt(t)/2).v,precision);
    return h*rt(t*(t+4))*(t-2)/(3*(t+1))+4*sq(h)/3-Real(11)*t/18+lp(t)/18+8*t/(9*(t+1));
}
struct State{int sign;unsigned weighted;std::array<std::array<int,4>,8> routes;std::vector<std::array<int,2>> U,F;};
struct Graph{std::string id,sha;std::vector<State> states;};
#include "mixed_single_weighted_e2_data.hpp"
Real determinant(const std::array<std::array<Real,4>,4>&matrix,unsigned size){
    std::array<unsigned,4> perm={0,1,2,3};Real total;
    do{unsigned inversions=0;Real term(1);for(unsigned i=0;i<size;++i){term=term*matrix[i][perm[i]];for(unsigned j=i+1;j<size;++j)inversions+=perm[i]>perm[j];}
       total=total+(inversions%2?-term:term);
    }while(std::next_permutation(perm.begin(),perm.begin()+size));return total;
}
Real polynomial(const std::vector<std::array<int,2>>&terms,const std::array<Real,8>&alpha){
    Real total;for(const auto&term:terms){Real value(term[1]);for(unsigned j=0;j<8;++j)if(term[0]&(1<<j))value=value*alpha[j];total=total+value;}return total;
}
std::array<Real,2> uf(const State&state,const std::array<Real,8>&alpha,bool positive){
    if(positive)return {polynomial(state.U,alpha),polynomial(state.F,alpha)};
    std::array<std::array<Real,4>,4>matrix;
    for(unsigned edge=0;edge<8;++edge)for(unsigned i=0;i<4;++i)for(unsigned j=0;j<4;++j)
        matrix[i][j]=matrix[i][j]+alpha[edge]*Real(state.routes[edge][i]*state.routes[edge][j]);
    return {determinant(matrix,3),determinant(matrix,4)};
}
Real density(const Graph&graph,const std::array<Real,8>&alpha,bool positive){
    Real mass;for(const auto&a:alpha)mass=mass+a;
    const unsigned weighted=graph.states.at(0).weighted;Real edge=alpha[weighted],t0=(mass-edge)/edge;
    pos(edge,"positive weighted edge parameter");pos(mass-edge,"remaining mass positive");
    Real a0=bubble_A(t0),c0=weighted_C(t0),total;
    for(const auto&state:graph.states){
        if(state.weighted!=weighted)throw std::runtime_error("fixed weighted edge required");
        auto values=uf(state,alpha,positive);pos(values[0],"positive U");pos(values[1],"positive F");
        Real q=values[1]/values[0],t1=t0+q/edge;
        Real J=edge*(t1*(bubble_A(t1)-a0)-(weighted_C(t1)-c0));
        total=total+Real(state.sign)*J/power3(values[0]);
    }
    return total;
}
Real value(const Graph&graph,bool positive,bool square,const std::array<Real,7>&p){
    std::array<Real,7>x;Real jacobian(1),total;
    for(unsigned j=0;j<7;++j){x[j]=square?sq(p[j]):p[j];if(square)jacobian=jacobian*2*p[j];}
    for(unsigned pivot=0;pivot<8;++pivot){std::array<Real,8>alpha;unsigned next=0;
        for(unsigned e=0;e<8;++e)alpha[e]=e==pivot?Real(1):x[next++];
        total=total+density(graph,alpha,positive);
    }
    total=total*jacobian;if(!arb_is_finite(total.v))throw std::runtime_error("finite complete E2 forest required");return total;
}
#include "mixed_schwinger_histogram_functions_v1.hpp"
Graph recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one weighted E2 recipe required");
    auto graph=load_weighted_e2(id);if(sha!=graph.sha)throw std::invalid_argument("exact original E2 identity required");return graph;}
}
int main(int argc,char**argv){try{using namespace weighted_e2;
    if(argc!=2)throw std::invalid_argument("one recipe required");auto graph=recipe(argv[1]);
    std::string mode;unsigned bits,count;if(!(std::cin>>mode>>bits>>count)||
      (mode!="hist-linear-original"&&mode!="hist-linear-positive"&&mode!="hist-square-original"&&mode!="hist-square-positive"&&mode!="map-linear"&&mode!="map-square")||bits<512||bits>2048||count<1||count>512)
        throw std::invalid_argument("bounded seven-dimensional weighted E2 request required");
    precision=bits;bool positive=mode.find("positive")!=std::string::npos,square=mode.find("square")!=std::string::npos;
    if(mode=="map-linear"||mode=="map-square"){
      std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"maps\":[";
      for(unsigned i=0;i<count;++i){std::array<Real,7>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict map cube");}
       auto m=schwinger_histogram_map(graph.id,square,p);if(i)std::cout<<',';std::cout<<"{\"point\":[";for(long j=0;j<7;++j){if(j)std::cout<<',';interval(m.point[j]);}std::cout<<"],\"jacobian\":";interval(m.jacobian);std::cout<<'}';}
      std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing map input");std::cout<<"]}\n";return 0;}
    std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
    for(unsigned i=0;i<count;++i){std::array<Real,7>p;for(auto&v:p){v=readq();if(!arb_is_positive(v.v)||!arb_lt(v.v,Real(1).v))throw std::invalid_argument("strict open seven-cube required");}
      auto result=schwinger_histogram_value(graph,positive,square,p);if(i)std::cout<<',';interval(result);}
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing weighted E2 request");std::cout<<"]}\n";return 0;
}catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 1;}}
