"""Diagnostic integration of the full host shape after analytic radial removal."""
import numpy as np
import primitive_child_triangle_host_radial_v1 as R
import primitive_host_grouped_stratified_vector_v1 as V


def norm(v):return (v[...,0]-v[...,1]/2)**2+3*v[...,1]**2/4
def dot(v,w):return (v[...,0]-v[...,1]/2)*(w[...,0]-w[...,1]/2)+3*v[...,1]*w[...,1]/4


def quadrature(order):
    assert order in (4,8,12)
    nodes,weights=np.polynomial.legendre.leggauss(order);nodes=(nodes+1)/2;weights=weights/2
    rows=[]
    for u,wu in zip(nodes,weights):
        for v,wv in zip(nodes,weights):rows.append(([u,(1-u)*v,(1-u)*(1-v)],wu*wv*(1-u)))
    assert abs(sum(w for _,w in rows)-.5)<1e-14 and all(w>0 for _,w in rows)
    return rows


def evaluate(points,binding,positive,proposal,order):
    p=np.asarray(points,dtype=np.longdouble);n=len(p);table=binding['groups']
    host=np.asarray(binding['host_routes'],dtype=np.longdouble)[:,1:];quad=quadrature(order)
    # The weights include the host triangle's area 1/2. This integrates
    # the full nonlinear radial kernel; no logarithm is evaluated at mean Q.
    hostdata=[]
    for coordinates,weight in quad:
        x=np.asarray(coordinates,dtype=np.longdouble);b=x@host;q=x[0]*x[1]+x[0]*x[2]+x[1]*x[2]
        hostdata.append((b,q,np.longdouble(weight)))
    total=np.zeros(n,dtype=np.longdouble)
    for st in proposal['strata']:
        a,x,t,tropical=V.parameters(p[:,:9],table,proposal,st)
        M=a.sum(axis=1);U=V.polynomial(table['U'],a);beta=np.zeros((n,9),dtype=np.longdouble)
        for g,edges in enumerate(table['groups']):
            if len(edges)==1:beta[:,edges[0]]=a[:,g]
            else:
                j=next(row['index'] for row in table['fractions'] if row['group']==g)
                beta[:,edges[0]]=a[:,g]*x[:,j];beta[:,edges[1]]=a[:,g]*(1-x[:,j])
        W=np.empty((n,3,3),dtype=np.longdouble)
        for i in range(3):
            for j in range(3):W[:,i,j]=V.polynomial(binding['child_variance_numerator'][i][j],beta)/(U*M)
        ell=W[:,0,0];assert (ell>0).all();ce=W[:,0,1:];energy=W[:,1,1]+W[:,2,2]-W[:,1,2]
        Ge=V.polynomial(positive['child_five_minor_Gram'],beta)/(U*M*M);assert (Ge>0).all()
        value=np.zeros(n,dtype=np.longdouble)
        for b,q,weight in hostdata:
            De=energy-2*dot(b,ce)+ell*norm(b);assert (De>0).all()
            # The prior control integrates to exact zero over this full
            # host domain; omitting it avoids a redundant cancellation.
            value+=weight*R.evaluate(np.full(n,q,dtype=np.longdouble),ell,De,Ge,np.zeros(n),False)
        factor=np.longdouble(str(float(st['probability']*table['normalization'])))
        total+=factor*(tropical/U)**3*value
    assert np.isfinite(total).all();return total
