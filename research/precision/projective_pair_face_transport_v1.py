"""Face-order proposals from exact paired numerator cancellations.

These are algebraic lower bounds on orders at a proper coordinate face.
They are not convergence or variance certificates. Existing forest densities
are unchanged; this module constructs a candidate importance distribution.
"""
from fractions import Fraction as Q
import math
import numpy as np
import exact_rational_sector_transport_v1 as H


def difference(left,right):
    values={int(m):int(c) for m,c in left}
    for mask,c in right:values[int(mask)]=values.get(int(mask),0)-int(c)
    return sorted([m,c] for m,c in values.items() if c)


def cross_difference(fb,uz,fz,ub):
    """Keep exponent multiplicities using union/intersection, not square-free masks."""
    values={}
    for sign,left,right in ((1,fb,uz),(-1,fz,ub)):
        for m,c in left:
            for n,d in right:
                key=(int(m)|int(n),int(m)&int(n))
                values[key]=values.get(key,0)+sign*int(c)*int(d)
    return sorted([u,v,c] for (u,v),c in values.items() if c)


def valuation(poly,face):
    if not poly:return None
    return min((face&int(u)).bit_count()+(face&int(v)).bit_count() for u,v,c in poly)


def all_valuations(poly,n):
    """Exact uint8 popcounts; batches bound temporary memory, with no floats."""
    if not poly:return None
    full=(1<<n)-1;assert n in (11,12)
    assert all(type(c) is int and c and 0<=v<=full and 0<u<=full and not(v&~u) for u,v,c in poly)
    union=np.array([u for u,v,c in poly],dtype=np.uint16)
    twice=np.array([v for u,v,c in poly],dtype=np.uint16)
    bits=np.array([i.bit_count() for i in range(1<<n)],dtype=np.uint8)
    out=np.zeros(full,dtype=np.int16)
    for start in range(1,full,64):
        masks=np.arange(start,min(start+64,full),dtype=np.uint16)
        orders=bits[union[:,None]&masks[None,:]]+bits[twice[:,None]&masks[None,:]]
        out[start:start+len(masks)]=orders.min(axis=0)
    return out


def polynomial(poly,alpha):
    return sum((Q(c)*math.prod(alpha[e] for e in range(len(alpha)) if u&(1<<e))*
                math.prod(alpha[e] for e in range(len(alpha)) if v&(1<<e)) for u,v,c in poly),Q())


def construct(binding,two_point):
    n=binding['original_edges'];assert n==(11 if two_point else 12)
    full=(1<<n)-1;states=binding['states'];vs=[];old=[];differences=[]
    for state in states:
        u=all_valuations([[m,0,int(c)] for m,c in state['U']],n)
        f=all_valuations([[m,0,int(c)] for m,c in state['F']],n)
        assert np.all(f>=u);vs.append((u,f));old.append(-3*u+(2 if two_point else 1)*(f-u))
    candidates=[];used=[]
    for b,z in binding['pairs']:
        assert states[b]['sign']==-states[z]['sign'];used.extend((b,z))
        ub,fb=vs[b];uz,fz=vs[z]
        du=[[m,0,c] for m,c in difference(states[b]['U'],states[z]['U'])]
        cross=cross_difference(states[b]['F'],states[z]['U'],states[z]['F'],states[b]['U'])
        vu=all_valuations(du,n);vc=all_valuations(cross,n)
        rz=fz-uz
        if vu is not None:
            inverse_difference=vu+2*np.minimum(ub,uz)-3*ub-3*uz
            candidates.append(inverse_difference+(2 if two_point else 1)*rz)
        if vc is not None:
            delta=vc-ub-uz;assert np.all(delta>=0)
            if two_point:
                candidates.append(-3*ub+delta+rz)
                candidates.append(-3*ub+2*delta)
            else:candidates.append(-3*ub+delta)
        differences.append({'pair':[b,z],'U_difference':du,'FU_cross_difference':cross})
    for index in binding.get('unpaired',[]):
        used.append(index);u,f=vs[index];candidates.append(-3*u+(2 if two_point else 1)*(f-u))
    assert sorted(used)==list(range(len(states))) and candidates
    cardinality=np.array([i.bit_count() for i in range(full)],dtype=np.int16)
    prior=np.min(np.stack(old),axis=0)+cardinality
    paired=np.min(np.stack(candidates),axis=0)+cardinality
    assert np.all(paired[1:]>=prior[1:]),'exact regrouping cannot lower the conservative face order'
    gamma={i:int(paired[i]) for i in range(1,full)}
    before={i:int(prior[i]) for i in range(1,full)}
    summary={'id':binding['id'],'two_point':two_point,'proper_faces':full-1,
        'improved_faces':sum(gamma[i]>before[i] for i in gamma),
        'critical_to_positive_faces':sum(before[i]<=0<gamma[i] for i in gamma),
        'remaining_nonpositive_faces':sum(gamma[i]<=0 for i in gamma),
        'minimum_old_gamma':min(before.values()),'minimum_new_gamma':min(gamma.values()),
        'U_difference_monomials':sum(len(d['U_difference']) for d in differences),
        'FU_cross_difference_monomials':sum(len(d['FU_cross_difference']) for d in differences),
        'density_changed':False,'absolute_integrability_certificate':False,'whole_integral_acceptance':False}
    return {'gamma':gamma,'old_gamma':before,'differences':differences,'summary':summary}


def transport(binding,two_point,power=2):
    result=construct(binding,two_point)
    # Cap positive proposal shapes at two, independently of numerical results.
    # For n<=12 and critical powers<=4 the crude integer bound
    # n!*4^(n-1)*2^(n-1) is below 2^62. This preserves the native table format.
    proposal={s:min(g,2) for s,g in result['gamma'].items()}
    assert math.factorial(binding['original_edges'])*8**(binding['original_edges']-1)<2**62
    result['proposal_positive_shape_cap']=2
    return H.build(binding['original_edges'],proposal,power),result
