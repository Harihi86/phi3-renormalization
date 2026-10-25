// The test driver substitutes the path of the already bound compiled source.
#define main phi3lab_shared_integrator_main
#include "__COMPILED_SOURCE__"
#undef main
int main(){try{using namespace weighted_e3;long bits,count;
 if(!(std::cin>>bits>>count)||bits<512||bits>2048||count<1||count>512)throw std::invalid_argument("bounded scalar request required");
 std::cout<<"{\"mode\":\"scalar-A\",\"bits\":"<<bits<<",\"count\":"<<count<<",\"values\":[";
 for(long i=0;i<count;++i){Rat t=read_rational();if(t<Rat(0))throw std::invalid_argument("nonnegative scalar argument required");Real value;bool done=false;std::string last;
  for(long working=bits;;working=std::min(working*2,8192L)){
   precision=working;
   try{value=shared_square_A(to_real(t));done=arb_is_finite(value.v)&&(arb_is_zero(value.v)||arb_rel_accuracy_bits(value.v)>=80);}
   catch(const std::runtime_error&e){last=e.what();}
   if(done)break;if(working==8192)throw std::runtime_error("scalar precision exhausted: "+last);
  }
  if(i)std::cout<<',';interval(value);
 }
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing scalar input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
