"""Floating diagnostics only: complete positive rank-two triangle update."""
from itertools import combinations,permutations
import numpy as np
import primitive_host_grouped_stratified_vector_v1 as V


def norm(v):return (v[...,0]-v[...,1]/2)**2+3*v[...,1]**2/4


def evaluate(points,binding,proposal,control=True,radial_reflection=2,simplex_orbit=1):
    p=np.asarray(points,dtype=np.longdouble);n=len(p);table=binding['groups']
    routes=np.asarray(binding['host_routes'],dtype=np.longdouble);loop=routes[:,:3];external=routes[:,3:]
    pairs=np.asarray(binding['triangle_pair_differences'],dtype=np.longdouble);L=pairs[:2,:3];E=pairs[:2,3:]
    cross=np.cross(L[0],L[1]);cofactors=[]
    for e,f in combinations(range(9),2):
        v=np.cross(loop[e],loop[f])
        if np.any(v):cofactors.append((e,f,v))
    root=np.sqrt(1-p[:,9]);x0=p[:,9]/(1+root)
    simplex=np.stack((x0,root*p[:,10],root*(1-p[:,10])),axis=1)
    orders=[(0,1,2)] if simplex_orbit==1 else list(permutations(range(3)))
    assert simplex_orbit in (1,6) and radial_reflection in (1,2)
    total=np.zeros(n,dtype=np.longdouble)
    for st in proposal['strata']:
        a,x,t,tropical=V.parameters(p[:,:9],table,proposal,st)
        assert ((t>0)&(t<1)).all()
        M=a.sum(axis=1);U=V.polynomial(table['U'],a);alpha=np.zeros((n,9),dtype=np.longdouble)
        center=np.zeros_like(alpha)
        for g,edges in enumerate(table['groups']):
            if len(edges)==1:alpha[:,edges[0]]=a[:,g];center[:,edges[0]]=a[:,g]
            else:
                j=next(row['index'] for row in table['fractions'] if row['group']==g)
                alpha[:,edges[0]]=a[:,g]*x[:,j];alpha[:,edges[1]]=a[:,g]*(1-x[:,j])
                center[:,edges[0]]=center[:,edges[1]]=a[:,g]/2
        adj=np.zeros((n,3,3),dtype=np.longdouble)
        for e,f,v in cofactors:adj+=alpha[:,e,None,None]*alpha[:,f,None,None]*np.outer(v,v)[None,:,:]
        inverse=adj/U[:,None,None]
        B=np.einsum('ne,ei,ej->nij',alpha,loop,external);Bc=np.einsum('ne,ei,ej->nij',center,loop,external)
        residual=E[None,:,:]-np.einsum('ai,nij,njk->nak',L,inverse,B)
        centered=E[None,:,:]-np.einsum('ai,nij,njk->nak',L,inverse,Bc)
        def covariance(v):
            return sum((alpha[:,e]*alpha[:,f]*np.dot(c,v)**2 for e,f,c in cofactors),np.zeros(n,dtype=np.longdouble))/U
        C00,C11,Cdiff=covariance(L[0]),covariance(L[1]),covariance(L[0]-L[1]);C01=(C00+C11-Cdiff)/2
        Cdet=np.einsum('ne,e->n',alpha,(loop@cross)**2)/U
        q0=V.polynomial(binding['host_F'],alpha)/U
        meanq=np.zeros(n,dtype=np.longdouble)
        if control:
            for mask,co in binding['host_F']:
                monomial=np.full(n,co,dtype=np.longdouble)
                for g,edges in enumerate(table['groups']):
                    count=sum(bool(mask&(1<<e)) for e in edges)
                    if len(edges)==2:
                        if count==1:monomial*=a[:,g]/2
                        elif count==2:monomial*=a[:,g]**2/6
                    elif count:monomial*=a[:,g]
                meanq+=monomial
            meanq/=U
        def adjquad(v,vector=True):
            if np.any(L[0]):
                ratio=C01/C00
                return C00*norm(v[:,1,:]-ratio[:,None]*v[:,0,:])+Cdet/C00*norm(v[:,0,:]) if vector else C00*(v[:,1]-ratio*v[:,0])**2+Cdet/C00*v[:,0]**2
            return C11*norm(v[:,0,:]) if vector else C11*v[:,0]**2
        adjpoint=adjquad(residual);adjcenter=adjquad(centered) if control else None
        zgroups=[]
        if control:
            for row in table['fractions']:
                l=loop[row['edges'][0]];z=np.einsum('ai,nij,j->na',L,inverse,l)
                zgroups.append((row['group'],z,adjquad(z,False)))
        value=np.zeros(n,dtype=np.longdouble)
        for order in orders:
            child=simplex[:,order];c01=child[:,0]*child[:,1];c02=child[:,0]*child[:,2];c12=child[:,1]*child[:,2]
            detK=np.prod(child,axis=1);traceK=c01*C00+c02*C11+c12*Cdiff
            primary=c01*norm(residual[:,0,:])+c02*norm(residual[:,1,:])+c12*norm(residual[:,0,:]-residual[:,1,:])
            if control:
                centerprimary=c01*norm(centered[:,0,:])+c02*norm(centered[:,1,:])+c12*norm(centered[:,0,:]-centered[:,1,:])
            for aux in ([t] if radial_reflection==1 else [t,1-t]):
                radial=M*aux/(1-aux);mass=M/(1-aux)
                determinant_delta=radial*traceK+radial**2*detK*Cdet;ratio=1+determinant_delta;inv3=ratio**-3
                one_minus=determinant_delta*(3+3*determinant_delta+determinant_delta**2)*inv3
                delta=radial*(primary+radial*detK*adjpoint)/ratio
                result=one_minus*np.log1p(q0/mass)-inv3*np.log1p(delta/(mass+q0))
                if control:
                    meandelta=radial*(centerprimary+radial*detK*adjcenter)/ratio
                    for g,z,adjz in zgroups:
                        primaryz=c01*z[:,0]**2+c02*z[:,1]**2+c12*(z[:,0]-z[:,1])**2
                        meandelta+=a[:,g]**2*radial*(primaryz+radial*detK*adjz)/(12*ratio)
                    result+=(one_minus*(meanq-q0)-inv3*(meandelta-delta))/mass
                value+=result/(aux*(1-aux))
        # Uniform simplex map has exact constant Jacobian 1/2.
        total+=np.longdouble(str(float(st['probability']*table['normalization'])))*(tropical/U)**3*value/(2*len(orders)*radial_reflection)
    assert np.isfinite(total).all();return total
