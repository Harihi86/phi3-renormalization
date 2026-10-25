// Exact simplex boundary recurrences, included after general scalar functions.
struct GeneralTriangleJet{Real value;std::array<Real,3>gradient;std::array<std::array<Real,3>,3>hessian;};
std::array<Real,4> general_edge_inverse(const Real&t){
    pos(t,"strict general edge invariant");Real as;arb_asinh(as.v,(rt(t)/2).v,precision);Real j=4*as/rt(t*(t+4));
    return {j,j/2,j/2+(j-1)/t,j/2+3*(j-1)/(2*t)};
}
GeneralTriangleJet general_triangle_jet(const Real&a,const Real&b,const Real&c){
    pos(a,"general jet a");pos(b,"general jet b");pos(c,"general jet c");
    Real d=(a+b-c)/2,det=a*b-sq(d);pos(det,"strict general jet Gram");
    std::array<std::array<Real,2>,2>inv;inv[0]={b/det,-d/det};inv[1]={-d/det,a/det};
    std::array<Real,2>u={(b*a-d*b)/(2*det),(a*b-d*a)/(2*det)};
    Real height=1+(a*u[0]+b*u[1])/2,J=inverse_triangle_general(a,b,c);
    Real Ea=edge_log_full(a),Eb=edge_log_full(b),Ec=edge_log_full(c);
    Real L=(u[0]*Eb+u[1]*Ea+(1-u[0]-u[1])*Ec)/2-Real(1)/2+height*J;
    std::array<Real,2>first;
    for(int i=0;i<2;++i)first[i]=u[i]*J-(inv[i][0]*(Ec-Eb)+inv[i][1]*(Ec-Ea))/2;
    std::array<std::array<Real,2>,2>B,K,second;B[0]={Ec/2,(Ec-Ea)/2};B[1]={(Ec-Eb)/2,Ec/2};
    std::array<Real,2>av={a,b};
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)K[i][j]=first[i]*av[j]+(i==j?L:Real(0))-B[i][j];
    for(int i=0;i<2;++i)for(int j=0;j<2;++j)for(int k=0;k<2;++k)second[i][j]=second[i][j]+K[i][k]*inv[k][j]/2;
    std::array<std::array<Real,5>,5>M,N;M[0][0]=J;M[1][0]=first[0];M[0][1]=first[1];M[2][0]=second[0][0];M[1][1]=second[0][1];M[0][2]=second[1][1];
    if(!arb_overlaps(second[0][1].v,second[1][0].v))throw std::runtime_error("inconsistent symmetric reciprocal moment");
    auto ea=general_edge_inverse(a),eb=general_edge_inverse(b),ec=general_edge_inverse(c);
    N[0][0]=(u[0]*eb[0]+u[1]*ea[0]+(1-u[0]-u[1])*ec[0])/(2*height);
    bool seen[5][5]={{false}};seen[0][0]=true;
    const int choose[4][4]={{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}};
    for(int degree=0;degree<4;++degree)for(int m=0;m<=degree;++m){int n=degree-m;Real bx,by;
        for(int k=0;k<=m;++k)bx=bx+Real((k%2?-1:1)*choose[m][k])*ec[n+k];
        for(int k=0;k<=n;++k)by=by+Real((k%2?-1:1)*choose[n][k])*ec[m+k];
        if(m==0)bx=bx-eb[n];if(n==0)by=by-ea[m];
        std::array<Real,2>v={bx-(m?Real(m)*M[m-1][n]:Real(0)),by-(n?Real(n)*M[m][n-1]:Real(0))};
        for(int i=0;i<2;++i){int mm=m+(i==0),nn=n+(i==1);Real next=u[i]*N[m][n]+(inv[i][0]*v[0]+inv[i][1]*v[1])/2;
            if(seen[mm][nn]){if(!arb_overlaps(N[mm][nn].v,next.v))throw std::runtime_error("inconsistent simplex recurrence paths");}
            else{N[mm][nn]=next;seen[mm][nn]=true;}}
    }
    const int terms[3][3][3]={{{1,0,1},{2,0,-1},{1,1,-1}},{{0,1,1},{1,1,-1},{0,2,-1}},{{1,1,1},{0,0,0},{0,0,0}}};
    GeneralTriangleJet out;out.value=L;
    for(int i=0;i<3;++i)for(int k=0;k<3;++k){const auto&left=terms[i][k];out.gradient[i]=out.gradient[i]+Real(left[2])*M[left[0]][left[1]];
        for(int j=0;j<3;++j)for(int l=0;l<3;++l){const auto&right=terms[j][l];out.hessian[i][j]=out.hessian[i][j]-Real(left[2]*right[2])*N[left[0]+right[0]][left[1]+right[1]];}}
    return out;
}
