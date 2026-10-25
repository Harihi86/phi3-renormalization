// One local T2 bubble and two complete T0 triangles, retaining the host loop.
struct BubbleTwoTrianglesGraph{std::string id,sha;std::array<std::array<int,3>,3>host;std::array<std::array<std::array<int,3>,3>,2>triangles;std::array<int,3>bubble;};
#include "one_bubble_two_triangles_data.hpp"
Real bubble_two_triangles_value(const BubbleTwoTrianglesGraph&g,bool positive,bool control,const std::array<Real,3>&p){
    Real s=p[0],z=p[1],theta=p[2],t=sq(s/(1-s)),A=1+t,co,si;
    arb_sin_cos_pi(si.v,co.v,theta.v,precision);Real rad=rt(t*z),sqrt3=rt(Real(3)),T0=-edge_log_full(t)/2,total;
    auto bj=mixed_bubble_jets(t);Real J=bj[0],P0=J*sq(T0);
    for(long sign:{-1L,1L}){
        Real qx=Real(sign)*rad*co,qy=Real(sign)*rad*si,q1=qx,q2=(-qx+sqrt3*qy)/2;
        auto norm=[&](const std::array<int,3>&r){
            if(positive)return sq(Real(r[0])*qx+Real(r[1])-Real(r[2])/2)+sq(Real(r[0])*qy+sqrt3*Real(r[2])/2)+Real(r[0]*r[0])*t*(1-z);
            return Real(r[0]*r[0])*t+2*Real(r[0])*(Real(r[1])*q1+Real(r[2])*q2)+Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        Real lr;for(const auto&r:g.host){Real d=1+norm(r);pos(d,"positive one-bubble host");lr=lr-lp((d-A)/A);}
        Real db=mixed_bubble_value(norm(g.bubble))-J;
        std::array<Real,2>dt;for(int k=0;k<2;++k)dt[k]=-triangle_log_full(norm(g.triangles[k][0]),norm(g.triangles[k][1]))-T0;
        Real difference=db*sq(T0)+(J+db)*(T0*(dt[0]+dt[1])+dt[0]*dt[1]);
        total=total+(em1(lr)*P0+real_exp(lr)*difference)/2;
    }
    Real common=2/(s*(1-s))*power3(t/A),answer=common*(1-z)*total;
    if(control){
        auto en=[](const std::array<int,3>&r){return Real(r[1]*r[1]+r[2]*r[2]-r[1]*r[2]);};
        auto lv=[](const std::array<int,3>&r){return std::array<Real,2>{Real(2*r[0]*r[1]),Real(2*r[0]*r[2])};};
        OneHostJet product;product.value=1;
        for(const auto&r:g.host){OneHostJet one;one.value=1;auto l=lv(r);one.constant=-en(r)/A;
            for(int i=0;i<2;++i){one.linear[i]=-l[i]/A;for(int j=0;j<2;++j)one.square[i][j]=l[i]*l[j]/sq(A);}product=multiply_one_host(product,one);}
        OneHostJet bubble;bubble.value=J;auto l=lv(g.bubble);bubble.constant=bj[1]*en(g.bubble);
        for(int i=0;i<2;++i){bubble.linear[i]=bj[1]*l[i];for(int j=0;j<2;++j)bubble.square[i][j]=bj[2]*l[i]*l[j]/2;}product=multiply_one_host(product,bubble);
        Real as;arb_asinh(as.v,(rt(t)/2).v,precision);Real j0=4*as/rt(t*(t+4)),j2=2*(1+j0)/(t+4),ell=-(1-j0)/(2*t),Cp=ell/2,Cf=-(j0-3*(1-j0)/t)/12,I2=(1-2*j0+j2)/sq(t);
        for(const auto&rows:g.triangles){
            auto a=lv(rows[0]),b=lv(rows[1]);for(int i=0;i<2;++i){a[i]=a[i]/2;b[i]=b[i]/2;}
            OneHostJet triangle;triangle.value=T0;triangle.constant=Cp*(en(rows[0])+en(rows[1]))+Cf*en(rows[2]);
            for(int i=0;i<2;++i){triangle.linear[i]=ell*(a[i]+b[i]);for(int j=0;j<2;++j)triangle.square[i][j]=I2*(a[i]*a[j]+b[i]*b[j]+(a[i]*b[j]+b[i]*a[j])/2)/3;}
            product=multiply_one_host(product,triangle);
        }
        std::array<Real,2>q={rad*co,rad*(-co+sqrt3*si)/2};Real point=product.constant;
        for(int i=0;i<2;++i)for(int j=0;j<2;++j)point=point+product.square[i][j]*q[i]*q[j];
        Real mean=product.constant+t*(product.square[0][0]+product.square[1][1]-(product.square[0][1]+product.square[1][0])/2)/6;
        answer=answer+common*(mean/2-(1-z)*point);
    }
    if(!arb_is_finite(answer.v))throw std::runtime_error("finite bubble-two-triangle forest required");return answer;
}
BubbleTwoTrianglesGraph bubble_two_triangles_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one original recipe required");auto g=load_bubble_two_triangles(id);if(sha!=g.sha)throw std::invalid_argument("original one-bubble identity required");return g;}
