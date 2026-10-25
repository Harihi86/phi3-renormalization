"""Independent unreduced720-order Gaussian reference, with both ordered weights."""
from fractions import Fraction as Q
from itertools import permutations
import mpmath as mp
import shared_e3_full_sector_reference_v1 as G
import nested_analytic_marginals_v1 as A
import bubble_positive_series_reference_v2 as B


def coefficients(point,binding,table):
    assert len(point)==6 and all(type(v) is Q and 0<v<1 for v in point)
    assert binding['weighted_kind']=='double'
    positive={};derivative={};weights=binding['weighted_edges'];count=0
    def put(target,t,value):target[t]=target.get(t,Q())+value
    auxiliary=point[5];y=1-auxiliary**2;s=(1-y*y)/4;measure=auxiliary*y**4/6
    for order in permutations(range(6)):
        remaining=set(range(6));remaining.remove(order[0]);alpha=[Q()]*6;alpha[order[0]]=Q(1);scale=Q(1);jac=Q(1)
        for j,u in enumerate(point[:5],1):
            mask=sum(2**edge for edge in remaining);shape=table['shapes'][mask]
            if shape<0:ratio=u**(-shape);density=(-shape)*u**(-shape-1)
            else:ratio=shape*u/(1-u+shape*u);density=Q(shape)/(1-u+shape*u)**2
            scale*=ratio;alpha[order[j]]=scale;jac*=density*ratio**(5-j);remaining.remove(order[j])
        assert not remaining;mass=sum(alpha,Q());d,e=alpha[weights[0]],alpha[weights[1]]
        for state in binding['paired_states']:
            assert state['weighted_edges']==weights
            U,z=G.schur(state,alpha);common=Q(state['sign'])*jac*measure/U**3
            # Expand the original c*J - h*e*dJ/dM + k*J_shift algebra directly.
            for m,prefactor in ((mass,(s-2)/(s-1)**2),(mass+e*(1/s-1),1/(s*(s-1)**2))):
                left=m/d-1;right=(m+z)/d-1
                put(positive,left,common*prefactor);put(positive,right,-common*prefactor)
            left=mass/d-1;right=(mass+z)/d-1
            put(derivative,left,-common*e*d/((1-s)*mass**2))
            put(derivative,right,common*e*d/((1-s)*(mass+z)**2))
        count+=1
    positive={t:c for t,c in positive.items() if c};derivative={t:c for t,c in derivative.items() if c}
    assert count==720 and sum(positive.values(),Q())==0 and all(t>0 for t in set(positive)|set(derivative))
    return positive,derivative


def evaluate_coefficients(pair,digits=640):
    ctx=mp.mp.clone();ctx.dps=digits;number=lambda q:ctx.mpf(q.numerator)/q.denominator
    a,b=pair;values=[]
    for t in sorted(set(a)|set(b)):
        x=number(t)
        if t in a:values.append(number(a[t])*A.bubble_primitive(ctx,x))
        if t in b:values.append(number(b[t])*B.hyper_bubble(ctx,x))
    return ctx,ctx.fsum(values)


def evaluate(point,binding,table,digits=640):
    return evaluate_coefficients(coefficients(point,binding,table),digits)
