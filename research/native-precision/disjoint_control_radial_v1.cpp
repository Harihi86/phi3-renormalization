// Exact angular/simplex control integral; finite ACB enclosure and explicit tails.
#define main frozen_mixed_triangle_main
#include "mixed-triangle-angular-native-v2/integrate.cpp"
#undef main
namespace disjoint_control {
using namespace mixed_triangle;
struct Params {unsigned graph;unsigned long calls=0,analytic=0,rejected=0,unsupported=0;};
Z pow(const Z& z,unsigned n){Z r;acb_pow_ui(r.v,z.v,n,precision);return r;}
int density(acb_ptr out,const acb_t x,void*opaque,slong order,slong prec){
    Params&p=*static_cast<Params*>(opaque);++p.calls;if(order)++p.analytic;
    if(order>1){++p.unsupported;for(slong j=0;j<order;++j)acb_indeterminate(out+j);return 0;}
    Z t;acb_exp(t.v,x,prec);
    // A positive real part throughout the input ball certifies all branches.
    if(!arb_is_positive(acb_realref(t.v))){++p.rejected;acb_indeterminate(out);return 0;}
    Z r,s,lg,tp=t+4,arg;acb_sqrt_analytic(r.v,t.v,order!=0,prec);
    acb_sqrt_analytic(s.v,tp.v,order!=0,prec);arg=(r+s)/2;
    acb_log_analytic(lg.v,arg.v,order!=0,prec);
    Z I0=4*lg/(r*s),J=tp*I0/2-2,I1=(1-I0)/t,J00=2*(1+I0)/tp;
    Z I2=(1-2*I0+J00)/square(t),A=1+t,B=tp*J/6-t/9;
    Z Ss=p.graph==32?Z(1):Z(Real(5)/2),Dot=p.graph==32?Z(Real(3)/2):Z(3);
    Z F=(5*t/(3*square(A))-2/A)*J/2+I0/12-t*Ss*I2/18-t*Dot*I1/(6*A);
    if(p.graph==42)F=F+I1/4;
    Z value=-pow(t,3)*B*F/(2*pow(A,4));acb_set(out,value.v);return 0;
}
}
int main(){try{using namespace disjoint_control;
    std::string mode;unsigned bits,cap;disjoint_control::Params p;
    if(!(std::cin>>mode>>bits>>cap>>p.graph)||bits<384||bits>2048||cap<16||cap>200000||
       (p.graph!=32&&p.graph!=42)||(mode!="point"&&mode!="piece"&&mode!="complex"&&mode!="tails"))
        throw std::invalid_argument("bounded control-integral request required");
    precision=bits;Real lo=readq(),hi=readq();
    if(!arb_ge(lo.v,Real(-24).v)||!arb_le(hi.v,Real(96).v)||!arb_le(lo.v,hi.v))
        throw std::invalid_argument("ordered natural-log domain [-24,96] required");
    Z a(lo),b(hi),result;
    if(mode=="complex"){
        Real rl=readq(),rh=readq(),il=readq(),ih=readq();
        if(!arb_le(rl.v,rh.v)||!arb_le(il.v,ih.v))throw std::invalid_argument("ordered complex box required");
        arb_union(acb_realref(a.v),rl.v,rh.v,bits);arb_union(acb_imagref(a.v),il.v,ih.v,bits);
    }
    std::string extra;if(std::cin>>extra)throw std::invalid_argument("trailing request input");
    if(mode=="tails"){
        Real low=arb_cw::exp(Real(-120))/150;
        Real high=Real(5)*arb_cw::exp(Real(-96))*Real(96*96+4*96+5)/3;
        std::cout<<"{\"low\":";interval(low);std::cout<<",\"high\":";interval(high);std::cout<<"}\n";return 0;
    }
    if(mode!="piece"){
        disjoint_control::density(result.v,a.v,&p,mode=="complex"?1:0,bits);
        std::cout<<"{\"finite\":"<<(acb_is_finite(result.v)?"true":"false")<<",\"rejected\":"<<p.rejected;
        if(acb_is_finite(result.v))value_json(result);std::cout<<"}\n";return 0;
    }
    if(!arb_cw::less(lo,hi))throw std::invalid_argument("positive piece width required");
    mag_t tol;mag_init(tol);mag_set_ui_2exp_si(tol,1,-136);
    acb_calc_integrate_opt_t options;acb_calc_integrate_opt_init(options);
    options->deg_limit=140;options->eval_limit=cap;options->depth_limit=1024;options->use_heap=1;
    int status=acb_calc_integrate(result.v,disjoint_control::density,&p,a.v,b.v,120,tol,options,bits);
    bool ok=status==ARB_CALC_SUCCESS&&acb_is_finite(result.v)&&arb_contains_zero(acb_imagref(result.v))&&!p.unsupported;
    std::cout<<"{\"calc_exit\":"<<status<<",\"bits\":"<<bits<<",\"graph\":"<<p.graph<<",\"calls\":"<<p.calls
        <<",\"analytic_calls\":"<<p.analytic<<",\"rejected_domains\":"<<p.rejected;
    if(acb_is_finite(result.v))value_json(result);
    std::cout<<",\"success\":"<<(ok?"true":"false")<<"}\n";mag_clear(tol);return ok?0:2;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
