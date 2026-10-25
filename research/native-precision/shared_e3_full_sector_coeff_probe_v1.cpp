// Includes a bound copy with exactly its top-level entry renamed by the test driver.
#include "__COMPILED_SOURCE__"
int main(int argc,char**argv){try{using namespace weighted_e3;
 if(argc!=2)throw std::invalid_argument("one bound recipe required");auto graph=recipe(argv[1]);std::string mode;long count;
 if(!(std::cin>>mode>>count)||count<1||count>16||(mode!="sector2-original"&&mode!="sector2-positive"&&mode!="sector4-original"&&mode!="sector4-positive"))throw std::invalid_argument("bounded exact coefficient request required");
 bool fourth=mode.find("sector4")!=std::string::npos,positive=mode.find("positive")!=std::string::npos;
 std::cout<<"{\"id\":\""<<graph.id<<"\",\"mode\":\""<<mode<<"\",\"count\":"<<count<<",\"coefficients\":[";
 for(long i=0;i<count;++i){std::array<Rat,5>point;for(auto&v:point){v=read_rational();if(!(Rat(0)<v)||!(v<Rat(1)))throw std::invalid_argument("five open exact coordinates required");}
  if(i)std::cout<<',';std::cout<<'[';bool first=true;
  for(const auto&entry:prepare_complete_shared(graph,point,positive,fourth)){
   if(!first)std::cout<<',';first=false;std::cout<<"[\""<<entry.first.str()<<"\",\""<<entry.second.str()<<"\"]";
  }
  std::cout<<']';
 }
 std::string extra;if(std::cin>>extra)throw std::invalid_argument("no trailing coefficient input");std::cout<<"]}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
