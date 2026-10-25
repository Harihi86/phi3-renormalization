"""Original Gaussian coefficients plus direct x integration; no closed primitive."""
from fractions import Fraction as Q


def routes_from_raw(raw):
    gd=raw['gaussian_full_descriptor'];out=[]
    for word in ([],[0],[1]):
        term=next(t for t in gd['terms'] if t['subgraphs_inner_first']==word)
        assert term['jet_subgraphs']==[] and term['jet_orders']==[]
        rows=[]
        for e in (0,1,3,4,6):
            old=term['scaled_routing'][e];assert all(not c['powers'] for c in old)
            values=[Q(c['coefficient']) for c in old]
            assert all(q+values[2]==0 for q in values[3:])
            rows.append(values[:3])
        out.append(rows)
    return out


def coefficients(a,b,y,routes):
    result=[]
    for rows in routes:
        values=[];U=None
        for x in (Q(),Q(1,2),Q(1)):
            alpha=(a*x,a*(1-x),b*y,b*(1-y),Q(1))
            m=[[sum((aa*row[i]*row[j] for aa,row in zip(alpha,rows)),Q()) for j in range(3)] for i in range(3)]
            uu=m[0][0]*m[1][1]-m[0][1]*m[1][0]
            ff=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
            assert uu>0
            if U is None:U=uu
            assert U==uu;values.append(ff/uu)
        c,mid,end=values;aa=2*(2*mid-c-end);bb=end-c+aa;result.append((U,aa,bb,c))
    return result


def evaluate(mp,gid,points,tp,control,routes):
    s,r,v,y,nu=map(Q,points);y/=2;nu/=2;t=(s/(1-s))**tp
    mq=lambda q:mp.mpf(q.numerator)/q.denominator
    mt=mq(t);density=mp.mpf(0);h=lambda z:(1+z)*mp.log1p(z)-z
    for a,b,jac in ((r*v,r,2*r),(r,1/v,2/v**2),(1/r,1/(r*v),2/(r**3*v**2))):
        rows=coefficients(a,b,y,routes);mass=mq(1+a+b)
        def integrand(x):
            terms=[]
            for u,aa,bb,c in rows:
                variance=-mq(aa)*x*x+mq(bb)*x+mq(c)
                terms.append(h(mt*variance/mass)/mq(u)**3)
            return terms[0]-terms[1]-terms[2]
        core=mp.quad(integrand,[0,mp.mpf('0.125'),mp.mpf('.5'),mp.mpf('.875'),1])
        density+=mq(jac*a*b*(1+a+b))*core
    if control:
        z0=mt*mq(y*(1-y)/(1+v))
        def limit(x):
            z1=mt*mq(v/(1+v))*x*(1-x)
            return h(z0+z1)-h(z0)-h(z1)
        density+=mq(2*(1+v)/v*(Q(1,control+1)-(1-r)**control))*mp.quad(limit,[0,mp.mpf('.5'),1])
    A=1+mt;AA=2+mt
    if gid=='L3E2G7':
        angular=mp.hyp2f1(1,mp.mpf('.5'),3,4*mt/(AA*AA))/AA
        outer=angular-(1/A-1/A**2+mp.mpf(2)/3*mt/A**3)
    else:
        assert gid=='L3E3G38'
        norm=mq(1-nu*(1-nu));angular=mp.hyp2f1(1,mp.mpf('1.5'),3,4*mt*norm/(AA*AA))/AA**2
        outer=angular-1/A**2
    measure=mq(t*t/(2*(1+t)**2)*tp*t/(s*(1-s)))
    return measure*outer*density
