"""Floating diagnostics for the complete primitive-child radial forest kernel."""
from itertools import permutations
from fractions import Fraction as Q
import numpy as np
import primitive_child_triangle_host_radial_v1 as R
import primitive_host_grouped_stratified_vector_v1 as V


def norm(v):return (v[...,0]-v[...,1]/2)**2+3*v[...,1]**2/4
def dot(v,w):return (v[...,0]-v[...,1]/2)*(w[...,0]-w[...,1]/2)+3*v[...,1]*w[...,1]/4


def proposal(table):
    gamma={int(k):v for k,v in table['gamma'].items()};nu={int(k):v for k,v in table['nu'].items()};strata=[]
    for order in permutations(range(6)):
        mask=63;mass=Q(1);exponents=[];leading=[]
        for e in order[:-1]:
            mask^=1<<e;mass/=gamma[mask];exponents.append(gamma[mask]);leading.append(nu[mask])
        strata.append({'order':list(order),'exponents':exponents,'leading':leading,'probability':mass/table['normalization']})
    assert len(strata)==720 and sum((st['probability'] for st in strata),Q())==1
    return {'kind':'all720','strata':strata}


def literal_parameters(p,table,proposal,st):
    ratios=p[:,1:6]**(1/np.asarray(st['exponents'],dtype=np.longdouble));a=np.ones((len(p),6),dtype=np.longdouble)
    a[:,st['order'][1:]]=np.cumprod(ratios,axis=1)
    return a,p[:,6:9],p[:,0],np.prod(ratios**np.asarray(st['leading']),axis=1)


def evaluate(points,binding,positive,proposal,control=True,radial_transport=True,simplex_orbit=6):
    p=np.asarray(points,dtype=np.longdouble);n=len(p);table=binding['groups']
    host=np.asarray(binding['host_routes'],dtype=np.longdouble)[:,1:];center=host.mean(axis=0)
    root=np.sqrt(1-p[:,9]);simplex=np.stack((p[:,9]/(1+root),root*p[:,10],root*(1-p[:,10])),axis=1)
    assert simplex_orbit in (1,6);orders=[(0,1,2)] if simplex_orbit==1 else list(permutations(range(3)))
    total=np.zeros(n,dtype=np.longdouble)
    for st in proposal['strata']:
        a,x,t,tropical=literal_parameters(p[:,:9],table,proposal,st);assert ((t>0)&(t<1)).all()
        M=a.sum(axis=1);U=V.polynomial(table['U'],a);beta=np.zeros((n,9),dtype=np.longdouble)
        for g,edges in enumerate(table['groups']):
            if len(edges)==1:beta[:,edges[0]]=a[:,g]
            else:
                j=next(row['index'] for row in table['fractions'] if row['group']==g)
                beta[:,edges[0]]=a[:,g]*x[:,j];beta[:,edges[1]]=a[:,g]*(1-x[:,j])
        W=np.empty((n,3,3),dtype=np.longdouble)
        for i in range(3):
            for j in range(3):W[:,i,j]=V.polynomial(binding['child_variance_numerator'][i][j],beta)/(U*M)
        ell=W[:,0,0];assert (ell>0).all()
        ce=W[:,0,1:];energy=W[:,1,1]+W[:,2,2]-W[:,1,2]
        Ge=V.polynomial(positive['child_five_minor_Gram'],beta)/(U*M*M);assert (Ge>0).all()
        meanDe=energy-2*dot(center,ce)+ell*norm(center)+ell/12
        value=np.zeros(n,dtype=np.longdouble)
        for order in orders:
            hostx=simplex[:,order];b=hostx@host;q0=hostx[:,0]*hostx[:,1]+hostx[:,0]*hostx[:,2]+hostx[:,1]*hostx[:,2]
            De=energy-2*dot(b,ce)+ell*norm(b);assert (De>0).all()
            value+=R.evaluate(q0,ell,De,Ge,meanDe,control)
        factor=np.longdouble(str(float(st['probability']*table['normalization'])))
        total+=factor*(tropical/U)**3*value/(2*len(orders))
    assert np.isfinite(total).all();return total
