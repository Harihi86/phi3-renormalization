"""Exact conditional subtraction, with an independent original-polynomial path."""
from fractions import Fraction as Q
import sympy as S
import two_loop_vertex_integrability_v1 as B
import two_loop_vertex_scale_integrated_v1 as M
F=B.F


def parameters(point,perm):
    assert len(point)==7 and all(type(x) in (int,Q) and 0<x<1 for x in point)
    s,v,f0,f1,f2,r,t=point;den=1+s+s*v;w=[None]*3
    for i,j in enumerate(perm):w[j]=(1/den,s/den,s*v/den)[i]
    fs=(f0,f1,f2);psi=w[0]*w[1]+w[0]*w[2]+w[1]*w[2];kap=w[0]*w[1]*w[2]/psi
    h=[ww*f*(1-f) for ww,f in zip(w,fs)];x,y=f0-f1,f0-f2
    a=h[0]+h[1]+kap*x*x;c=h[0]+kap*x*y
    determinant=h[0]*h[1]+h[0]*h[2]+h[1]*h[2]+kap*(h[0]*(y-x)**2+h[1]*y*y+h[2]*x*x)
    host=[r,(1-r)*t,(1-r)*(1-t)];rho=c/a
    q=host[0]*host[1]+host[0]*host[2]+host[1]*host[2]
    n=a*F.norm(rho-host[0]-host[2],rho-host[0])
    return {'Q':q,'u':a,'N':n,'d':determinant/a,'average_N':a*(rho*rho-rho+Q(5,12)),
        'average_Q':Q(1,4),'host_jacobian':1-r,'weight':v/(1+v+s*v)**3,
        'alpha':[w[i]*fs[i] for i in range(3)]+[w[i]*(1-fs[i]) for i in range(3)]+host,
        'projective_weight':w[0]*w[1]*w[2]*s/den**3,'child_U':psi}


def polynomial_parameters(raw,point,perm):
    base=parameters(point,perm);states=raw['positive_full_descriptor']['states']
    bare=next(s for s in states if s['word']==[]);ct=next(s for s in states if s['word']==[0])
    assert bare['sign']==1 and ct['sign']==-1
    def extract(alpha):
        u=B.coefficients(bare['U'],alpha,(2,3));f=B.coefficients(bare['F'],alpha,(2,3,4))
        uc=B.coefficients(ct['U'],alpha,(2,));fc=B.coefficients(ct['F'],alpha,(2,))
        assert u[2]==uc[2] and f[2]==fc[2]
        q=f[2]/u[2];a=u[3]/u[2];d=f[4]/u[3];n=f[3]/u[2]-a*q-d
        return {'Q':q,'u':a,'N':n,'d':d,'U0':u[2]}
    result=extract(base['alpha']);centres=[]
    for zero in range(3):
        alpha=base['alpha'][:6]+[Q() if j==zero else Q(1,2) for j in range(3)]
        p=extract(alpha);assert p['u']==result['u'] and p['d']==result['d'];centres.append(p)
    result.update({'average_N':sum((r['N'] for r in centres),Q())/3,
        'average_Q':sum((r['Q'] for r in centres),Q())/3,
        'host_jacobian':base['alpha'][7]+base['alpha'][8],
        'weight':base['projective_weight']/result.pop('U0')**3})
    assert result=={k:base[k] for k in result}
    return result


def evaluate(data,mp):
    def m(v):return mp.mpf(v.numerator)/v.denominator if type(v) is Q else mp.mpf(v)
    q,u,n,d,an,aq,hj,weight=[m(data[k]) for k in ('Q','u','N','d','average_N','average_Q','host_jacobian','weight')]
    H=[-mp.log(u)/(1-u)]
    for k in range(1,4):H.append((H[-1]-mp.mpf(1)/k)/(1-u))
    a=u*sum(H[:3]);first=q*a-n*H[3]-d*H[2];average=aq*a-an*H[3]-d*H[2]
    return weight*(hj*(M.scale_integral(q,u,n,d,mp)-first)+average/2)


def orbit(point,perm,mp):
    return sum((evaluate(parameters(tuple(point[:2])+tuple(point[2+i] for i in fp)+tuple(point[5:]),perm),mp) for fp in F.PERMUTATIONS),mp.mpf(0))/6


def symbolic():
    u,t,sigma=S.symbols('u t sigma',positive=True)
    H=[]
    for k in range(1,5):
        primitive=-S.log(1-(1-u)*t)/(1-u)**k-sum(t**j/(j*(1-u)**(k-j)) for j in range(1,k))
        assert S.simplify(S.diff(primitive,t)-t**(k-1)/(1-(1-u)*t))==0
        H.append(S.simplify(primitive.subs(t,1)))
    for k in range(1,4):assert S.simplify(H[k]-(H[k-1]-S.Rational(1,k))/(1-u))==0
    R=1+u*sigma
    assert S.cancel((1-R**-3)/sigma-u*(R**-1+R**-2+R**-3))==0
    h0,h2,rho=S.symbols('h0 h2 rho');h1=1-h0-h2
    q=h0*h1+h0*h2+h1*h2;n=F.norm(rho-h0-h2,rho-h0)
    def average(poly):
        pp=S.Poly(S.expand(poly),h0,h2)
        return S.expand(sum(co*2*S.factorial(i)*S.factorial(j)/S.factorial(i+j+2) for (i,j),co in pp.terms()))
    assert average(q)==S.Rational(1,4)
    assert S.simplify(average(n)-(rho*rho-rho+S.Rational(5,12)))==0
    for i in range(3):
        for j in range(3-i):
            poly=h0**i*h2**j
            edge_average=sum(poly.subs({h0:a,h2:b}) for a,b in ((S.Rational(1,2),S.Rational(1,2)),(0,S.Rational(1,2)),(S.Rational(1,2),0)))/3
            assert average(poly)==edge_average
    return {'H_antiderivatives_exact':4,'H_recurrences_exact':3,'sigma_pair_identity_exact':True,
        'host_Q_average':'1/4','host_N_over_u_average':'rho^2-rho+5/12',
        'independent_midpoint_rule_exact_for_total_degree_leq_2':6,
        'identity':'weight*[(1-r)*(J-J1)+average(J1)/2]',
        'full_J_retained_no_momentum_truncation':True,
        'f_permutation_orbit':{'count':6,'each_coefficient':'1/6','unit_absolute_jacobian':True}}
