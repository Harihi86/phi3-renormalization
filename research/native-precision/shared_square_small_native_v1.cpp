#include <cmath>
#include <chrono>
#include "shared_full_sector_frozen.cpp"
namespace weighted_e3 {
#include "shared_square_small_native_v1.hpp"
}
int main(){try{using namespace weighted_e3;
 std::string mode;long bits,count;
 if(!(std::cin>>mode>>bits>>count)||(mode!="series"&&mode!="closed")||bits<512||bits>2048||count<1||count>128)throw std::invalid_argument("bounded scalar comparison required");
 std::cout<<"{\"mode\":\""<<mode<<"\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(long i=0;i<count;++i){Rat input=read_rational();if(input<Rat(0)||Rat("1/16")<input)throw std::invalid_argument("zero-to-one-sixteenth argument required");
  precision=bits;auto tick=std::chrono::steady_clock::now();SmallSquareResult result;
  if(mode=="series")result=shared_square_small(to_real(input));
  else{for(long work=bits;;work=std::min(2*work,8192L)){
   precision=work;Real value=shared_square_A(to_real(input));
   if(arb_is_finite(value.v)&&(arb_is_zero(value.v)||arb_rel_accuracy_bits(value.v)>=80)){result={value,Real(0),0,false,work};break;}
   if(work==8192)throw std::runtime_error("closed scalar relative precision exhausted");
  }}
  if(!arb_is_zero(result.value.v)&&arb_rel_accuracy_bits(result.value.v)<80)throw std::runtime_error("minimum scalar relative bits required");
  double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
  if(i)std::cout<<',';std::cout<<"{\"input\":\""<<input.str()<<"\",\"series\":"<<(result.series?"true":"false")<<",\"terms\":"<<result.terms<<",\"working_bits\":"<<result.working_bits<<",\"seconds\":"<<elapsed<<",\"interval\":";interval(result.value);
  std::cout<<",\"tail\":";interval(result.tail);std::cout<<'}';
 }
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing scalar input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
