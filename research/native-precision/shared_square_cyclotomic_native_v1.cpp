// The original complete-sector executable is compiled unchanged as a reference.
#include <cmath>
#include <chrono>
#define main frozen_complete_sector_main
#include "shared_full_sector_frozen.cpp"
#undef main
namespace weighted_e3 {
#include "shared_square_cyclotomic_native_v1.hpp"
}
int main(){try{using namespace weighted_e3;
 std::string mode;long bits,count;
 if(!(std::cin>>mode>>bits>>count)||(mode!="series"&&mode!="closed")||bits<512||bits>2048||count<1||count>128)throw std::invalid_argument("bounded scalar comparison required");
 precision=bits;std::cout<<"{\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(long i=0;i<count;++i){Rat input=read_rational();if(!(Rat(0)<input)||!(input<Rat(1)))throw std::invalid_argument("strict positive unit interval required");
  Real w=to_real(input);auto tick=std::chrono::steady_clock::now();CyclicResult result;
  if(mode=="series")result=shared_square_cyclic(w);else result={shared_square_P(w,false),Real(0),0,false};
  double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
  if(i)std::cout<<',';std::cout<<"{\"input\":\""<<input.str()<<"\",\"series\":"<<(result.series?"true":"false")<<",\"terms\":"<<result.terms<<",\"seconds\":"<<elapsed<<",\"interval\":";interval(result.value);
  std::cout<<",\"tail\":";interval(result.tail);std::cout<<'}';
 }
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing scalar input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
