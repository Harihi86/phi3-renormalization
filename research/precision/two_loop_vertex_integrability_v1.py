"""Exact positive identities and original-polynomial scale extraction for G35."""
from fractions import Fraction as Q
import sympy as S
import two_loop_vertex_rank_one_v1 as F


def symbolic():
    w=S.symbols('w0:3',positive=True);f=S.symbols('f0:3',positive=True)
    psi=w[0]*w[1]+w[0]*w[2]+w[1]*w[2];kap=S.prod(w)/psi
    h=[w[i]*f[i]*(1-f[i]) for i in range(3)]
    x,y=f[0]-f[1],f[0]-f[2]
    a=h[0]+h[1]+kap*x*x;b=h[0]+h[2]+kap*y*y;c=h[0]+kap*x*y
    positive=[]
    for i in range(3):
        j,k=[n for n in range(3) if n!=i]
        # w_i-kappa = w_i^2 (w_j+w_k)/Psi.
        assert S.cancel(w[i]-kap-w[i]**2*(w[j]+w[k])/psi)==0
        positive.append((w[i]-kap)*f[i]*(1-f[i])+kap*(f[i]*(1-f[j])*(1-f[k])+(1-f[i])*f[j]*f[k]))
    assert all(S.cancel(left-right)==0 for left,right in zip((c,a-c,b-c),positive))
    sigma,ps,aa,qq,nn,dd=S.symbols('sigma psi a Q N d',positive=True)
    u=ps*sigma**2*(1+sigma*aa)
    ff=ps*sigma**2*(qq+sigma*(aa*qq+nn+dd)+sigma**2*aa*dd)
    u0=S.expand(u).coeff(sigma,2);u1=S.expand(u).coeff(sigma,3)
    f0=S.expand(ff).coeff(sigma,2);f1=S.expand(ff).coeff(sigma,3);f2=S.expand(ff).coeff(sigma,4)
    extracted=(f0/u0,u1/u0,f1/u0-u1*f0/u0**2-f2/u1,f2/u1)
    assert all(S.cancel(l-r)==0 for l,r in zip(extracted,(qq,aa,nn,dd)))
    return {'nonnegative_forms':['c','a-c','b-c'],'port_ratio':'0 <= c/a <= 1',
        'bounds':['0<Q<=1','0<u<=1','0<=N/u<=3','0<=d<=1',
                  'u>=s*v*f0*(1-f0)/3'],
        'absolute_sigma_majorant':'7/3+2*log(4)-log(s)-log(v)-log(f0)-log(1-f0)',
        'untransformed_remaining_measure_bound':'v on each of six unit cubes',
        'original_scale_coefficients_exact':True,'rigorous_quadrature_error_bound':False}


def coefficients(poly,alpha,allowed):
    ans={degree:Q() for degree in allowed}
    for row in poly:
        mask,slot,channel=row[:3];assert slot==0
        degree=(mask&63).bit_count();assert degree in allowed
        value=Q(row[3]) if len(row)==4 else Q(1)
        for i in range(9):
            if mask&(1<<i):value*=alpha[i]
        ans[degree]+=value
    return ans


def extract(raw,point,perm):
    chart=F.chart(tuple(point)+(Q(1,2),),perm)
    alpha=[v/chart['sigma'] if i<6 else v for i,v in enumerate(chart['alphas'])]
    states=raw['positive_full_descriptor']['states']
    bare=next(s for s in states if s['word']==[]);ct=next(s for s in states if s['word']==[0])
    assert bare['sign']==1 and ct['sign']==-1
    u=coefficients(bare['U'],alpha,(2,3));ff=coefficients(bare['F'],alpha,(2,3,4))
    uc=coefficients(ct['U'],alpha,(2,));fc=coefficients(ct['F'],alpha,(2,))
    assert u[2]==uc[2] and ff[2]==fc[2]
    q=ff[2]/u[2];a=u[3]/u[2];d=ff[4]/u[3];n=ff[3]/u[2]-a*q-d
    result={'Q':q,'u':a,'N':n,'d':d}
    assert result=={k:chart[k] for k in result}
    assert 0<q<=1 and 0<a<=1 and 0<=n/a<=3 and 0<=d<=1
    s,v,f0=point[:3];assert a>=s*v*f0*(1-f0)/3
    return result
