"""Exact forest differences after integrating the complete weighted bubble.

The density is sum sign*J(F/U; M,d)/U**3 on eight effective parameters.
Only the proposal changes. Face orders allow logarithms and are not a joint
multiscale convergence, variance, or whole-integral error certificate.
"""
from fractions import Fraction as Q
import math
import numpy as np
import projective_pair_face_transport_v1 as P
import projective_pair_face_transport_v2 as V
import exact_rational_sector_transport_v1 as H

N=8
FULL=(1<<N)-1


def valuations(poly):
    if not poly:return None
    assert all(type(c) is int and c and 0<=v<=FULL and 0<u<=FULL and not(v&~u) for u,v,c in poly)
    masks=np.arange(FULL,dtype=np.uint16)
    bits=np.array([m.bit_count() for m in range(FULL+1)],dtype=np.uint8)
    out=np.full(FULL,1000,dtype=np.int16)
    for u,v,c in poly:out=np.minimum(out,bits[masks&u]+bits[masks&v])
    return out


def validate(binding,option):
    states=binding['paired_states'];used=[]
    for b,z in option['pairs']:
        assert 0<=b<len(states) and 0<=z<len(states)
        x,y=set(states[b]['word']),set(states[z]['word'])
        assert x<y and len(y-x)==1 and states[b]['sign']==-states[z]['sign']
        used.extend((b,z))
    used.extend(option['unpaired'])
    assert sorted(used)==list(range(len(states)))


def construct(binding):
    assert binding['effective_edge_count']==N and binding['effective_loops']==3
    states=binding['paired_states'];weighted=binding['weighted_edge']
    assert all(s['weighted_edge']==weighted for s in states)
    options=V.matchings(states);assert options
    orders=[];single=[];old=[]
    # b=M-d can vanish on just the face containing all unweighted edges.
    nb=np.array([int((s|(1<<weighted))==FULL) for s in range(FULL)],dtype=np.int16)
    for state in states:
        u=valuations([[int(m),0,int(c)] for m,c in state['U']])
        f=valuations([[int(m),0,int(c)] for m,c in state['F']]);assert np.all(f>=u)
        h=f-u;orders.append((u,h));old.append(-3*u+2*h)
        single.append(-3*u+2*h+2*np.minimum(nb,h))
    edges=sorted(set(tuple(p) for option in options for p in option['pairs']))
    bounds={};numerators=[]
    for b,z in edges:
        ub,hb=orders[b];uz,hz=orders[z]
        du=[[m,0,c] for m,c in P.difference(states[b]['U'],states[z]['U'])]
        cross=P.cross_difference(states[b]['F'],states[z]['U'],states[z]['F'],states[b]['U'])
        vu,vc=valuations(du),valuations(cross);pieces=[]
        if vu is not None:pieces.append(vu+2*np.minimum(ub,uz)-3*ub-3*uz+2*hz+2*np.minimum(nb,hz))
        if vc is not None:
            delta=vc-ub-uz;assert np.all(delta>=0)
            pieces.append(-3*ub+delta+hz+2*np.minimum(nb,hz))
            pieces.append(-3*ub+2*delta+2*np.minimum(np.minimum(nb,hz),delta))
        bounds[(b,z)]=np.min(np.stack(pieces),axis=0) if pieces else np.full(FULL,1000,dtype=np.int16)
        numerators.append({'pair':[b,z],'U_difference':du,'FU_cross_difference':cross})
    sizes=np.array([m.bit_count() for m in range(FULL)],dtype=np.int16)
    original=np.min(np.stack(old),axis=0)+sizes
    best=np.full(FULL,-1000,dtype=np.int16);witness=np.full(FULL,-1,dtype=np.int16)
    for i,option in enumerate(options):
        validate(binding,option)
        pieces=[bounds[tuple(pair)] for pair in option['pairs']]+[single[j] for j in option['unpaired']]
        value=np.min(np.stack(pieces),axis=0)+sizes
        improve=value>best;best=np.maximum(best,value);witness[improve]=i
    assert np.all(best[1:]>=original[1:]) and np.all(best[1:]<1000)
    gamma={s:int(best[s]) for s in range(1,FULL)}
    return {'id':binding['id'],'gamma':gamma,'original_gamma':{s:int(original[s]) for s in gamma},
        'matchings':options,'witness':{s:int(witness[s]) for s in gamma},'numerators':numerators,
        'weighted_edge':weighted,'summary':{'id':binding['id'],'proper_faces':FULL-1,
        'complete_matchings':len(options),'original_forests':len(binding['states']),
        'bubble_grouped_states':len(states),'improved_faces':int(sum(gamma[s]>original[s] for s in gamma)),
        'critical_to_positive_faces':int(sum(original[s]<=0<gamma[s] for s in gamma)),
        'remaining_nonpositive_faces':sum(gamma[s]<=0 for s in gamma),
        'minimum_gamma':min(gamma.values()),'whole_graph_acceptance':False}}


def transport(binding,cap):
    assert cap in (2,3)
    proof=construct(binding)
    assert min(proof['gamma'].values())>0,'positive proposal requires every face witness'
    table=H.build(N,{s:min(g,cap) for s,g in proof['gamma'].items()},2)
    return table,proof


def polynomial(poly,alpha):
    return sum((Q(c)*math.prod(alpha[e] for e in range(N) if mask&(1<<e)) for mask,c in poly),Q())


def gaussian_uf(state,alpha):
    routes=state['routes']
    a=[[sum((v*r[i]*r[j] for v,r in zip(alpha,routes)),Q()) for j in range(4)] for i in range(4)]
    u=Q(1)
    for k in range(3):
        pivot=a[k][k];assert pivot>0;u*=pivot
        for i in range(k+1,4):
            for j in range(k+1,4):a[i][j]-=a[i][k]*a[k][j]/pivot
    return u,u*a[3][3]


def coefficients(binding,alpha,gaussian=False):
    """Combine all A/C coefficients exactly, including common lower endpoints."""
    assert len(alpha)==N and all(type(a) in (int,Q) and a>0 for a in alpha)
    d=alpha[binding['weighted_edge']];t0=(sum(alpha)-d)/d;values={}
    def add(kind,t,co):values[(kind,t)]=values.get((kind,t),Q())+co
    for state in binding['paired_states']:
        u,f=gaussian_uf(state,alpha) if gaussian else (polynomial(state['U'],alpha),polynomial(state['F'],alpha))
        assert u>0 and f>0
        t1=t0+f/(d*u);co=state['sign']*d/u**3
        add(0,t1,co*t1);add(1,t1,-co);add(0,t0,-co*t1);add(1,t0,co)
    return [[kind,t,co] for (kind,t),co in sorted(values.items()) if co]
