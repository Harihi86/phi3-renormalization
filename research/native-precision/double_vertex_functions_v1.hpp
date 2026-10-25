// The two parent branches are combined before any integration.
struct DoubleVertexGraph{std::string id,sha;std::array<bool,2>bubble;std::array<std::array<int,4>,3>outer,inner;std::array<std::array<std::array<int,4>,3>,2>child;};
#include "double_vertex_data.hpp"
Real dv_baseline(const std::array<int,4>&r,const Real&t,const Real&v,const Real&h){return Real(r[0]*r[0])*t+Real(r[1]*r[1])*v+2*Real(r[0]*r[1])*h;}
QuadraticJet dv_triangle_jet(const std::array<std::array<int,4>,3>&rows,const Real&t,const Real&v,const Real&h){
    QuadraticJet child;std::array<Real,3>base;for(int a=0;a<3;++a)base[a]=dv_baseline(rows[a],t,v,h);
    if(rows[2][0]==0&&rows[2][1]==0){
        if(rows[0][0]!=rows[1][0]||rows[0][1]!=rows[1][1]||!arb_equal(external_norm4(rows[2]).v,Real(1).v))throw std::invalid_argument("exact unit nested triangle required");
        Real x=base[0],as;pos(x,"positive unit triangle baseline");arb_asinh(as.v,(rt(x)/2).v,precision);
        Real j0=4*as/rt(x*(x+4)),j2=2*(1+j0)/(x+4),ell=-(1-j0)/(2*x),cp=ell/2,cf=-(j0-3*(1-j0)/x)/12,i2=(1-2*j0+j2)/sq(x);
        V4 a=linear_vector4(rows[0]),b=linear_vector4(rows[1]);child.value=-edge_log_full(x)/2;child.linear=scale4(sum4(a,b),ell);
        child.constant=cp*(external_norm4(rows[0])+external_norm4(rows[1]))+cf*external_norm4(rows[2]);
        for(int i=0;i<4;++i)for(int j=0;j<4;++j)child.square[i][j]=i2*(a[i]*a[j]+b[i]*b[j]+(a[i]*b[j]+b[i]*a[j])/2)/3;
    }else{
        auto jet=general_triangle_jet(base[0],base[1],base[2]);child.value=-jet.value;
        std::array<V4,3>vectors={linear_vector4(rows[0]),linear_vector4(rows[1]),linear_vector4(rows[2])};
        for(int a=0;a<3;++a){child.linear=sum4(child.linear,scale4(vectors[a],-2*jet.gradient[a]));child.constant=child.constant-jet.gradient[a]*external_norm4(rows[a]);
            for(int b=0;b<3;++b)for(int i=0;i<4;++i)for(int j=0;j<4;++j)child.square[i][j]=child.square[i][j]-2*jet.hessian[a][b]*vectors[a][i]*vectors[b][j];}
    }return child;
}
Real dv_bubble(const Real&x){pos(x,"positive shifted nested bubble");return ((x+4)*edge_log_full(x)/6-x/9)/(1+x);}
QuadraticJet dv_bubble_jet(const std::array<int,4>&row,const Real&t,const Real&v,const Real&h){
    Real x=dv_baseline(row,t,v,h),as;pos(x,"positive nested bubble baseline");arb_asinh(as.v,(rt(x)/2).v,precision);
    Real e=edge_log_full(x),j=4*as/rt(x*(x+4)),b=(x+4)*e/6-x/9,b1=e/6+(x+4)*(1-j)/(6*x)-Real(1)/9,b2=Real(1)/(6*x)-(1-j)/sq(x),d=1+x;
    Real first=b1/d-b/sq(d),second=b2/d-2*b1/sq(d)+2*b/power3(d);V4 linear=scale4(linear_vector4(row),Real(2));
    QuadraticJet result;result.value=b/d;result.linear=scale4(linear,first);result.constant=first*external_norm4(row);
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)result.square[i][j]=second*linear[i]*linear[j]/2;return result;
}
std::array<Real,2> double_vertex_chart(const DoubleVertexGraph&g,bool positive,const V4&q,const V4&r,int excluded){
    Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),w=scalar4(sum4(q,r),sum4(q,r));pos(t,"positive nested q");pos(v,"positive nested k");pos(w,"positive nested sum");
    std::array<Real,3>d={1+t,1+v,1+w};Real weight=power3(d[excluded])/(power3(d[0])+power3(d[1])+power3(d[2]));
    Real sqrt3=rt(Real(3));V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    auto baseline=[&](const std::array<int,4>&rr){return dv_baseline(rr,t,v,h);};
    std::array<QuadraticJet,2>O;
    Real K0=1/power3(1+v);
    for(int i=0;i<2;++i){O[i]=g.bubble[i]?dv_bubble_jet(g.child[i][0],t,v,h):dv_triangle_jet(g.child[i],t,v,h);
        K0=K0*(g.bubble[i]?dv_bubble(v):-edge_log_full(v)/2);}

    std::vector<std::array<int,4>>hosts(g.outer.begin(),g.outer.end());hosts.insert(hosts.end(),g.inner.begin(),g.inner.end());
    std::array<Real,2>den={Real(1),Real(1)},differences;std::array<QuadraticJet,2>jets;
    for(int branch=0;branch<2;++branch){
        int count=branch==0?6:3;jets[branch].value=1;
        for(int a=0;a<count;++a){const auto&row=hosts[a];Real dd=1+baseline(row);pos(dd,"positive nested host baseline");den[branch]=den[branch]*dd;
            QuadraticJet one;one.value=1;one.linear=scale4(linear_vector4(row),Real(-2)/dd);one.constant=-external_norm4(row)/dd;
            for(int i=0;i<4;++i)for(int j=0;j<4;++j)one.square[i][j]=one.linear[i]*one.linear[j];jets[branch]=multiply_quadratic(jets[branch],one);}
        if(branch==0)for(int i=0;i<2;++i)jets[branch]=multiply_quadratic(jets[branch],O[i]);
    }
    for(long sign:{-1L,1L}){
        auto norm=[&](const std::array<int,4>&rr){
            if(positive){V4 vec;for(int i=0;i<4;++i)vec[i]=Real(sign)*(Real(rr[0])*q[i]+Real(rr[1])*r[i]);vec[0]=vec[0]+Real(rr[2])-Real(rr[3])/2;vec[1]=vec[1]+sqrt3*Real(rr[3])/2;return scalar4(vec,vec);}
            return baseline(rr)+2*Real(sign)*scalar4(linear_vector4(rr),projection)+external_norm4(rr);};
        std::array<Real,2>os,od;
        for(int i=0;i<2;++i){os[i]=g.bubble[i]?dv_bubble(norm(g.child[i][0])):-triangle_log_general(norm(g.child[i][0]),norm(g.child[i][1]),norm(g.child[i][2]));od[i]=os[i]-O[i].value;}

        std::array<Real,2>lr;
        for(int a=0;a<6;++a){Real base=1+baseline(hosts[a]),shifted=1+norm(hosts[a]);pos(shifted,"positive shifted nested host");Real logratio=-lp((shifted-base)/base);lr[0]=lr[0]+logratio;if(a<3)lr[1]=lr[1]+logratio;}
        differences[0]=differences[0]+(em1(lr[0])*O[0].value*O[1].value+real_exp(lr[0])*(od[0]*O[1].value+O[0].value*od[1]+od[0]*od[1]))/2;
        differences[1]=differences[1]+em1(lr[1])/2;
    }
    std::array<Real,2>point,mean;
    for(int branch=0;branch<2;++branch){point[branch]=mean[branch]=jets[branch].constant;
        for(int i=0;i<4;++i)for(int j=0;j<4;++j){Real loop=(i/2==0&&j/2==0)?t:(i/2==1&&j/2==1)?v:h;Real external=i%2==j%2?Real(1):Real(-1)/2;
            point[branch]=point[branch]+jets[branch].square[i][j]*projection[i]*projection[j];mean[branch]=mean[branch]+jets[branch].square[i][j]*loop*external/6;}}
    Real chi=1/power3(1+sq(sq((t+v)/64)));
    Real complete=differences[0]/den[0]-K0*differences[1]/den[1],control=point[0]/den[0]-K0*point[1]/den[1],average=mean[0]/den[0]-K0*mean[1]/den[1];
    return {weight*(complete-chi*control),weight*chi*average};
}
DoubleVertexGraph double_vertex_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one nested graph recipe required");auto g=load_double_vertex(id);if(sha!=g.sha)throw std::invalid_argument("exact original nested identity required");return g;}
