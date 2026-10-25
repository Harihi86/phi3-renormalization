"""Original G5 effective Gaussian routes and independent hypergeometric bubble."""
from fractions import Fraction as Q


def routes_from_raw(raw):
    assert raw['id']=='L3E2G5';gd=raw['gaussian_full_descriptor'];result=[]
    for word in ([],[1],[2]):
        term=next(t for t in gd['terms'] if t['subgraphs_inner_first']==word)
        assert term['jet_subgraphs']==[]
        rows=[]
        for e in (1,7,2,5,0):
            coeff=term['scaled_routing'][e];assert all(not q['powers'] for q in coeff)
            row=[Q(q['coefficient']) for q in coeff];assert row[0]==0
            rows.append([row[1],-row[2],row[3]])
        result.append(rows)
    return result


def uf(a,b,x,y,route):
    alpha=(a*x,a*(1-x),b*y,b*(1-y),Q(1));m=[[sum((aa*row[i]*row[j] for aa,row in zip(alpha,route)),Q()) for j in range(3)] for i in range(3)]
    u=m[0][0]*m[1][1]-m[0][1]*m[1][0]
    f=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
    assert u>0;return u,f/u


def evaluate(mp,points,control,routes):
    r,v,x,y,u=map(Q,points);mq=lambda q:mp.mpf(q.numerator)/q.denominator
    bubble=lambda z:z*z/60*mp.hyp2f1(1,1,mp.mpf('3.5'),-z/4)
    total=mp.mpf(0)
    for aa,bb,jac in ((r*v,r,r),(r,1/v,1/v**2),(1/r,1/(r*v),1/(r**3*v**2))):
        for a,b in ((aa,bb),(bb,aa)):
            M=1+a+b;d=b*(1-y);point=Q();average=Q();primitives=[]
            for state,route in enumerate(routes):
                U,q=uf(a,b,x,y,route);sign=1 if state==0 else -1
                c=uf(a,b,Q(),y,route)[1];mid=uf(a,b,Q(1,2),y,route)[1];end=uf(a,b,Q(1),y,route)[1]
                A=2*(2*mid-c-end);B=end-c+A
                mean=A*A/5-A*B/2+(B*B-2*A*c)/3+B*c+c*c
                point+=sign*q*q/U**3;average+=sign*mean/U**3;primitives.append((U,q,sign))
            for aux in (u,1-u):
                value=mp.mpf(0)
                for U,q,sign in primitives:
                    mass=M+aux*q
                    value+=mq(sign*d*q*q*(1-aux)/(mass*mass*U**3))*bubble(mq((mass-d)/d))
                if control:value+=mq(d/(M*M)*(average/2-(1-aux)*point))*bubble(mq((M-d)/d))
                total+=mq(jac*a*b/2)*value
    return total
