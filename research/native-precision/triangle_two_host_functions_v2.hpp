// Two nested triangle loops integrated, retaining every parent forest branch.
struct TriangleTwoHostGraph{std::string id,sha;bool parent;std::array<std::array<int,4>,6>host;std::array<std::array<int,4>,3>inner,inner0;std::array<std::array<int,5>,3>tri,tri0;};
#include "triangle_two_host_data.hpp"
QuadraticJet tth_const(const Real&x){QuadraticJet a;a.value=x;return a;}
QuadraticJet tth_add(const QuadraticJet&a,const QuadraticJet&b){QuadraticJet c;c.value=a.value+b.value;c.constant=a.constant+b.constant;c.linear=sum4(a.linear,b.linear);for(int i=0;i<4;++i)for(int j=0;j<4;++j)c.square[i][j]=a.square[i][j]+b.square[i][j];return c;}
QuadraticJet tth_scale(const QuadraticJet&a,const Real&b){QuadraticJet c;c.value=a.value*b;c.constant=a.constant*b;c.linear=scale4(a.linear,b);for(int i=0;i<4;++i)for(int j=0;j<4;++j)c.square[i][j]=a.square[i][j]*b;return c;}
QuadraticJet tth_sub(const QuadraticJet&a,const QuadraticJet&b){return tth_add(a,tth_scale(b,Real(-1)));}
QuadraticJet tth_unary(const QuadraticJet&a,const Real&value,const Real&first,const Real&second){QuadraticJet c;c.value=value;c.linear=scale4(a.linear,first);c.constant=a.constant*first;for(int i=0;i<4;++i)for(int j=0;j<4;++j)c.square[i][j]=first*a.square[i][j]+second*a.linear[i]*a.linear[j]/2;return c;}
QuadraticJet tth_inv(const QuadraticJet&a){return tth_unary(a,1/a.value,-1/sq(a.value),2/power3(a.value));}
QuadraticJet tth_mul(const QuadraticJet&a,const QuadraticJet&b){return multiply_quadratic(a,b);}
QuadraticJet tth_div(const QuadraticJet&a,const QuadraticJet&b){return tth_mul(a,tth_inv(b));}
QuadraticJet tth_log(const QuadraticJet&a){return tth_unary(a,real_log(a.value),1/a.value,-1/sq(a.value));}
QuadraticJet tth_sqrt(const QuadraticJet&a){pos(a.value,"positive radial root discriminant");Real r=rt(a.value);return tth_unary(a,r,1/(2*r),-1/(4*a.value*r));}
std::array<Real,3> tth_master(const Real&t){
    pos(t,"positive radial master argument");Real pi;arb_const_pi(pi.v,precision);Real at1=sq(pi)/6-Real(5)/4,z=t-1,abs;
    arb_abs(abs.v,z.v);
    if(arb_lt(abs.v,(Real(1)/2).v)){
        mag_t radius,tail;mag_init(radius);mag_init(tail);arb_get_mag(radius,z.v);double upper=mag_get_d(radius);
        unsigned n=upper>0?unsigned(std::ceil((precision+40)/std::max(1.0,-std::log2(upper))))+3:3;n=std::max(3u,std::min(n,unsigned(precision+48)));
        std::vector<Real>powers(n+1);powers[0]=1;for(unsigned k=1;k<=n;++k)powers[k]=powers[k-1]*z;
        std::array<Real,3>out={at1,Real(0),Real(0)};
        for(unsigned k=1;k<=n;++k){Real c=Real((k%2)?2:-2)/(Real(long(k))*Real(long(k))*Real(long(k+1))*Real(long(k+2)));out[0]=out[0]+c*powers[k];out[1]=out[1]+Real(long(k))*c*powers[k-1];if(k>=2)out[2]=out[2]+Real(long(k))*Real(long(k-1))*c*powers[k-2];}
        for(unsigned j=0;j<3;++j){mag_pow_ui(tail,radius,n+1-j);mag_mul_2exp_si(tail,tail,2);arb_add_error_mag(out[j].v,tail);}mag_clear(radius);mag_clear(tail);return out;
    }
    Real u=1-t,log=real_log(t),li;Z dil=dilog(Z(1-t));if(!arb_contains_zero(acb_imagref(dil.v)))throw std::runtime_error("real radial dilogarithm required");arb_set(li.v,acb_realref(dil.v));
    Real value=sq(pi)/6-li+t*(2-3*t)*log/(2*sq(u))-t/(2*u);
    Real first=1/(2*u)-t/sq(u)-sq(t)*log/power3(u);
    Real second=-1/(2*sq(u))-t*(3+2*log)/power3(u)-3*sq(t)*log/sq(sq(u));return {value,first,second};
}
QuadraticJet tth_phi(const QuadraticJet&a){auto f=tth_master(a.value);return tth_unary(a,f[0],f[1],f[2]);}
QuadraticJet tth_radial(const std::array<QuadraticJet,3>&v,const Real&ell){
    auto one=tth_const(1),a=tth_add(one,v[0]),b=tth_add(tth_add(one,tth_scale(a,ell)),v[1]),c=tth_scale(tth_add(one,v[2]),ell);
    auto diff=tth_add(tth_sub(one,tth_scale(a,ell)),v[1]);auto disc=tth_add(tth_mul(diff,diff),tth_scale(tth_mul(a,tth_sub(v[1],v[2])),4*ell));
    auto root=tth_sqrt(disc),plus=tth_scale(tth_div(tth_add(b,root),a),Real(1)/2),minus=tth_scale(tth_div(c,tth_add(b,root)),2),L=tth_log(a);
    auto result=tth_add(tth_scale(L,real_log(ell)+Real(3)/2),tth_scale(tth_mul(L,L),Real(1)/2));
    result=tth_sub(tth_sub(result,tth_phi(tth_scale(plus,1/ell))),tth_phi(tth_scale(minus,1/ell)));
    Real pi;arb_const_pi(pi.v,precision);return tth_add(result,tth_const(tth_master(1/ell)[0]+sq(pi)/6-Real(5)/4));
}
Real tth_radial_value(const std::array<QuadraticJet,3>&v,const Real&ell,const V4&projection,long sign,bool shifted){
    std::array<QuadraticJet,3>args;for(int i=0;i<3;++i){Real x=v[i].value;if(shifted)x=x+Real(sign)*scalar4(v[i].linear,projection)+v[i].constant;args[i]=tth_const(x);}return tth_radial(args,ell).value;
}
V4 tth_row(const std::array<int,4>&r){return {Real(r[0]),Real(r[1]),Real(r[2]),Real(r[3])};}
QuadraticJet tth_norm(const V4&co,const Real&t,const Real&v,const Real&h){
    QuadraticJet a;a.value=sq(co[0])*t+sq(co[1])*v+2*co[0]*co[1]*h;a.constant=sq(co[2])+sq(co[3])-co[2]*co[3];a.linear={2*co[0]*co[2],2*co[0]*co[3],2*co[1]*co[2],2*co[1]*co[3]};return a;
}
std::array<QuadraticJet,3> tth_geometry(const std::array<std::array<int,4>,3>&inner,const std::array<std::array<int,5>,3>&tri,const std::array<Real,3>&a,const std::array<Real,3>&b,const Real&t,const Real&v,const Real&h,Real&ell){
    std::array<V4,3>rr,tt;V4 center;for(int i=0;i<3;++i){rr[i]=tth_row(inner[i]);center=sum4(center,scale4(rr[i],a[i]));for(int j=0;j<4;++j)tt[i][j]=Real(tri[i][j+1]);}
    std::array<Real,3>w={b[0]*b[1],b[0]*b[2],b[1]*b[2]};ell=0;int zero=-1;for(int i=0;i<3;++i){ell=ell+w[i]*Real(tri[i][0]*tri[i][0]);if(!tri[i][0])zero=i;}
    pos(ell,"positive inner triangle Gaussian coefficient");if(zero<0)throw std::runtime_error("one zero inner-loop difference required");
    std::array<QuadraticJet,3>out;int n=0;for(int i=0;i<3;++i)for(int j=i+1;j<3;++j){out[0]=tth_add(out[0],tth_scale(tth_norm(sum4(rr[i],rr[j],-1),t,v,h),a[i]*a[j]));++n;}
    for(int i=0;i<3;++i)out[1]=tth_add(out[1],tth_scale(tth_norm(sum4(tt[i],scale4(center,Real(-tri[i][0]))),t,v,h),w[i]));
    out[2]=tth_scale(tth_norm(tt[zero],t,v,h),b[0]*b[1]*b[2]/ell);return out;
}
std::array<Real,2> triangle_two_host_chart(const TriangleTwoHostGraph&g,bool positive,const V4&q,const V4&r,int excluded,const std::array<Real,3>&a,const std::array<Real,3>&b){
    Real t=scalar4(q,q),v=scalar4(r,r),h=scalar4(q,r),w=scalar4(sum4(q,r),sum4(q,r));pos(t,"positive triangle host q");pos(v,"positive triangle host k");pos(w,"positive triangle host sum");
    std::array<Real,3>d={1+t,1+v,1+w};Real weight=power3(d[excluded])/(power3(d[0])+power3(d[1])+power3(d[2])),sqrt3=rt(Real(3));
    V4 projection={q[0],(-q[0]+sqrt3*q[1])/2,r[0],(-r[0]+sqrt3*r[1])/2};
    Real ell,parent0;auto geometry=tth_geometry(g.inner,g.tri,a,b,t,v,h,ell);auto K=tth_radial(geometry,ell);Real K0=K.value;
    if(g.parent){Real ell0;auto zero=tth_geometry(g.inner0,g.tri0,a,b,t,v,h,ell0);parent0=tth_radial(zero,ell0).value/power3(1+v);}
    std::array<Real,2>den={Real(1),Real(1)},differences;std::array<QuadraticJet,2>jets;
    for(int branch=0;branch<(g.parent?2:1);++branch){int count=branch?3:6;jets[branch].value=1;
        for(int j=0;j<count;++j){auto norm=tth_norm(tth_row(g.host[j]),t,v,h);Real dd=1+norm.value;den[branch]=den[branch]*dd;auto prop=tth_scale(tth_inv(tth_add(tth_const(1),norm)),dd);jets[branch]=tth_mul(jets[branch],prop);}if(!branch)jets[branch]=tth_mul(jets[branch],K);}
    for(long sign:{-1L,1L}){
        Real full;
        if(positive){
            // Reconstruct actual four-component shifted vectors independently of dot expansion.
            auto raw_norm=[&](const V4&co){V4 x;for(int i=0;i<4;++i)x[i]=Real(sign)*(co[0]*q[i]+co[1]*r[i]);x[0]=x[0]+co[2]-co[3]/2;x[1]=x[1]+sqrt3*co[3]/2;return scalar4(x,x);};
            std::array<V4,3>rr,tt;V4 center;for(int i=0;i<3;++i){rr[i]=tth_row(g.inner[i]);center=sum4(center,scale4(rr[i],a[i]));for(int j=0;j<4;++j)tt[i][j]=Real(g.tri[i][j+1]);}
            std::array<Real,3>ww={b[0]*b[1],b[0]*b[2],b[1]*b[2]};std::array<QuadraticJet,3>gg;int zero=-1;
            for(int i=0;i<3;++i){for(int j=i+1;j<3;++j)gg[0].value=gg[0].value+a[i]*a[j]*raw_norm(sum4(rr[i],rr[j],-1));gg[1].value=gg[1].value+ww[i]*raw_norm(sum4(tt[i],scale4(center,Real(-g.tri[i][0]))));if(!g.tri[i][0])zero=i;}
            gg[2].value=b[0]*b[1]*b[2]/ell*raw_norm(tt[zero]);full=tth_radial(gg,ell).value;
        }else full=tth_radial_value(geometry,ell,projection,sign,true);
        std::array<Real,2>lr;for(int j=0;j<6;++j){auto co=tth_row(g.host[j]);auto norm=tth_norm(co,t,v,h);Real base=1+norm.value,shifted;
            if(positive){V4 x;for(int i=0;i<4;++i)x[i]=Real(sign)*(co[0]*q[i]+co[1]*r[i]);x[0]=x[0]+co[2]-co[3]/2;x[1]=x[1]+sqrt3*co[3]/2;shifted=1+scalar4(x,x);}else shifted=base+Real(sign)*scalar4(norm.linear,projection)+norm.constant;
            pos(shifted,"positive shifted triangle host");Real logratio=-lp((shifted-base)/base);lr[0]=lr[0]+logratio;if(j<3)lr[1]=lr[1]+logratio;}
        differences[0]=differences[0]+(em1(lr[0])*K0+real_exp(lr[0])*(full-K0))/2;if(g.parent)differences[1]=differences[1]+em1(lr[1])/2;
    }
    std::array<Real,2>point,mean;for(int branch=0;branch<(g.parent?2:1);++branch){point[branch]=mean[branch]=jets[branch].constant;
        for(int i=0;i<4;++i)for(int j=0;j<4;++j){Real loop=(i/2==0&&j/2==0)?t:(i/2==1&&j/2==1)?v:h;Real ext=i%2==j%2?Real(1):Real(-1)/2;point[branch]=point[branch]+jets[branch].square[i][j]*projection[i]*projection[j];mean[branch]=mean[branch]+jets[branch].square[i][j]*loop*ext/6;}}
    Real chi=1/power3(1+sq(sq((t+v)/64))),complete=differences[0]/den[0]-parent0*differences[1]/den[1],control=point[0]/den[0]-parent0*point[1]/den[1],average=mean[0]/den[0]-parent0*mean[1]/den[1];return {weight*(complete-chi*control),weight*chi*average};
}
TriangleTwoHostGraph triangle_two_host_recipe(const char*path){std::ifstream in(path);std::string id,sha,extra;if(!(in>>id>>sha)||(in>>extra))throw std::invalid_argument("one nested triangle recipe required");auto g=load_triangle_two_host(id);if(sha!=g.sha)throw std::invalid_argument("exact original nested triangle identity required");return g;}
