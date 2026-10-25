// Full nested two-loop child, with an optional three-loop parent counterterm.
struct BubbleTwoHostGraph{std::string id,sha;bool parent;std::array<std::array<int,4>,6>host;std::array<std::array<int,4>,3>child,parent0;};
#include "bubble_two_host_data.hpp"
Real bth_base(const std::array<int,4>&r,const Real&t,const Real&v,const Real&h){return Real(r[0]*r[0])*t+Real(r[1]*r[1])*v+2*Real(r[0]*r[1])*h;}
std::array<Real,2> bubble_two_host_chart(const BubbleTwoHostGraph&g,bool positive,const V4&q,const V4&r,int excluded,const Real&x,const Real&y){
    Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),w=scalar4(sum4(q,r),sum4(q,r));pos(t,"positive child host q");pos(v,"positive child host k");pos(w,"positive child host sum");
    std::array<Real,3>d={1+t,1+v,1+w},weights={x*(1-x)*y,x*(1-x)*(1-y),sq(1-x)*y*(1-y)};
    Real weight=power3(d[excluded])/(power3(d[0])+power3(d[1])+power3(d[2])),sqrt3=rt(Real(3));
    V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    auto baseline=[&](const std::array<int,4>&rr){return bth_base(rr,t,v,h);};
    Real V0,Cp,Vparent;V4 linear;
    for(int j=0;j<3;++j){V0=V0+weights[j]*baseline(g.child[j]);Cp=Cp+weights[j]*external_norm4(g.child[j]);linear=sum4(linear,scale4(linear_vector4(g.child[j]),2*weights[j]));
        if(g.parent)Vparent=Vparent+weights[j]*baseline(g.parent0[j]);}
    pos(V0,"positive child variance");Real arg=(1-x+V0)/x,lower=bubble_A((1-x)/x),at0=bubble_A(arg),B=bubble_value(arg),D=1+V0;
    Real K0=-(at0-lower),Kp=-x*B/sq(D),Kpp=2*x*B/power3(D)-bubble_derivative(arg)/sq(D),parent0;
    if(g.parent){pos(Vparent,"positive parent child variance");parent0=-(bubble_A((1-x+Vparent)/x)-lower)/power3(1+v);}
    QuadraticJet K;K.value=K0;K.constant=Kp*Cp;K.linear=scale4(linear,Kp);
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)K.square[i][j]=Kpp*linear[i]*linear[j]/2;
    std::array<Real,2>den={Real(1),Real(1)},differences;std::array<QuadraticJet,2>jets;
    for(int branch=0;branch<(g.parent?2:1);++branch){
        int count=branch==0?6:3;jets[branch].value=1;
        for(int a=0;a<count;++a){const auto&row=g.host[a];Real dd=1+baseline(row);pos(dd,"positive child host baseline");den[branch]=den[branch]*dd;
            QuadraticJet one;one.value=1;one.linear=scale4(linear_vector4(row),Real(-2)/dd);one.constant=-external_norm4(row)/dd;
            for(int i=0;i<4;++i)for(int j=0;j<4;++j)one.square[i][j]=one.linear[i]*one.linear[j];jets[branch]=multiply_quadratic(jets[branch],one);}
        if(branch==0)jets[branch]=multiply_quadratic(jets[branch],K);
    }
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){V4 vec;for(int i=0;i<4;++i)vec[i]=Real(sign)*(Real(rr[0])*q[i]+Real(rr[1])*r[i]);vec[0]=vec[0]+Real(rr[2])-Real(rr[3])/2;vec[1]=vec[1]+sqrt3*Real(rr[3])/2;return scalar4(vec,vec);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        Real V;for(int i=0;i<3;++i)V=V+weights[i]*norm(g.child[i]);pos(V,"positive shifted child variance");Real DK=-(bubble_A((1-x+V)/x)-at0);
        std::array<Real,2>lr;for(int a=0;a<6;++a){Real base=1+baseline(g.host[a]),shifted=1+norm(g.host[a]);pos(shifted,"positive shifted child host");Real logratio=-lp((shifted-base)/base);lr[0]=lr[0]+logratio;if(a<3)lr[1]=lr[1]+logratio;}
        differences[0]=differences[0]+(em1(lr[0])*K0+real_exp(lr[0])*DK)/2;
        if(g.parent)differences[1]=differences[1]+em1(lr[1])/2;
    }
    std::array<Real,2>point,mean;
    for(int branch=0;branch<(g.parent?2:1);++branch){point[branch]=mean[branch]=jets[branch].constant;
        for(int i=0;i<4;++i)for(int j=0;j<4;++j){Real loop=(i/2==0&&j/2==0)?t:(i/2==1&&j/2==1)?v:h;Real external=i%2==j%2?Real(1):Real(-1)/2;
            point[branch]=point[branch]+jets[branch].square[i][j]*projection[i]*projection[j];mean[branch]=mean[branch]+jets[branch].square[i][j]*loop*external/6;}}
    Real chi=1/power3(1+sq(sq((t+v)/64)));
    Real complete=differences[0]/den[0]-parent0*differences[1]/den[1],control=point[0]/den[0]-parent0*point[1]/den[1],average=mean[0]/den[0]-parent0*mean[1]/den[1];
    return {weight*(complete-chi*control),weight*chi*average};
}
BubbleTwoHostGraph bubble_two_host_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one nested child recipe required");auto g=load_bubble_two_host(id);if(sha!=g.sha)throw std::invalid_argument("exact original nested child identity required");return g;}
