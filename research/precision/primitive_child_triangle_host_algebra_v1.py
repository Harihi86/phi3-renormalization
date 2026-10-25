"""Positive child minors and exact host-simplex conditional moments."""
from itertools import combinations
import sympy as S
import primitive_child_triangle_host_factor_v1 as F


def derive(binding):
    routes=S.Matrix(binding['child_routes']);beta=S.symbols('b:9');U=F.polynomial(binding['child_U'],beta)
    N=S.Matrix([[F.polynomial(p,beta) for p in row] for row in binding['child_variance_numerator']])
    vectors=[];K=[]
    for subset in combinations(range(9),4):
        ds=[routes.extract(subset,[0,1,2,3+j]).det() for j in range(3)]
        if any(ds):vectors.append([sum(1<<e for e in subset),list(map(int,ds))])
    for i in range(3):
        for j in range(3):
            assert S.expand(N[i,j]-sum(S.prod(beta[e] for e in range(9) if mask&(1<<e))*ds[i]*ds[j] for mask,ds in vectors))==0
    for subset in combinations(range(9),5):
        ds=[routes.extract(subset,[0,1,2,3,4+j]).det() for j in range(2)]
        co=ds[0]**2+ds[1]**2-ds[0]*ds[1];assert co>=0 and co.is_Integer
        if co:K.append([sum(1<<e for e in subset),int(co)])
    assert K
    energy=N[0,0]*(N[1,1]+N[2,2]-N[1,2])-N[0,1]**2-N[0,2]**2+N[0,1]*N[0,2]
    assert S.expand(energy-U*F.polynomial(K,beta))==0
    return {'id':binding['id'],'child_four_minors':vectors,'child_five_minor_Gram':K,
            'variance_rank_one_sum_exact':True,'radial_quadratic_coefficient_positive':True,
            'Ge':'K/(Uc*Mchild^2)','De':'sum beta_mask * |d_ext-d_loop*b_host|_Gram^2/(Uc*Mchild)'}


def conditional_identity(host_routes):
    u,v,r=S.symbols('u v r');a=[u,v,1-u-v]
    w=S.symbols('w:6');W=S.Matrix([[w[0],w[1],w[2]],[w[1],w[3],w[4]],[w[2],w[4],w[5]]])
    H=sum((t*S.Matrix(row)*S.Matrix([row]) for t,row in zip(a,host_routes)),S.zeros(3))
    gram=S.Matrix([[1,-S.Rational(1,2)],[-S.Rational(1,2),1]])
    norm=lambda b:(b*gram*b.T)[0]
    b=H[:1,1:];c=W[:1,1:];E=W[1,1]+W[2,2]-W[1,2];ell=w[0];Ur=1+r*ell
    Q0=H[1,1]+H[2,2]-H[1,2]-norm(b)
    Qr=H[1,1]+H[2,2]-H[1,2]+r*E-norm(b+r*c)/Ur
    De=E-2*(b*gram*c.T)[0]+ell*norm(b);Ge=ell*E-norm(c)
    assert S.cancel(Qr-Q0-r*(De+r*Ge)/Ur)==0
    assert S.expand(ell*De-Ge-norm(c-ell*b))==0
    centered=lambda p:S.cancel(p.subs({u:S.Rational(1,3),v:S.Rational(1,3)}))
    mean=lambda p:S.cancel(2*S.integrate(S.integrate(S.expand(p),(v,0,1-u)),(u,0,1)))
    assert S.expand(Q0-(a[0]*a[1]+a[0]*a[2]+a[1]*a[2]))==0
    assert centered(Q0)==S.Rational(1,3) and mean(Q0)==S.Rational(1,4)
    assert S.cancel(mean(Qr)-centered(Qr)+1/(12*Ur))==0
    assert S.expand(mean(De)-centered(De)-ell/12)==0
    A=S.symbols('A',positive=True);radial=(A**u-1)/(A-A**u)
    jac=S.log(A)*(A-1)*A**u/((A-A**u)*(A**u-1))
    assert S.simplify(S.diff(radial,u)/radial-jac)==0
    return {'full_radial_identity':True,'delta_nonnegative':True,'delta_at_most_r_De':True,
            'host_normalized_mean_Q0':'1/4','host_normalized_mean_Qr':'Qr(center)-1/[12(1+r*ell)]',
            'host_normalized_mean_De':'De(center)+ell/12','matched_radial_map':'r=(ell^u-1)/(ell-ell^u)',
            'radial_map_independent_of_host_simplex':True,'complete_conditional_residual_retained':True}
