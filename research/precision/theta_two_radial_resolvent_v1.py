"""Apply the proven normalized Gauss resolvent bound to two radial directions.

No native code, complex region, external expansion, tail or allocation changes.
Native cell output retains its conservative historical bound. A separate
certificate encloses its actual central Gauss sum using the sharper exact
remainder; a changed order still requires a freshly evaluated central sum.
"""
from fractions import Fraction as Q
import theta_two_radial_majorant_audit_v5 as A
import gauss_cauchy_resolvent as RES


def bounds(box,radii):
    assert len(box)==4 and len(radii)==2
    xl,xr,yl,yr=map(lambda q:A.bounded(Q(q)),box)
    rx,ry=map(lambda q:A.bounded(Q(q)),radii)
    assert xl<xr and yl<yr and rx>(xr-xl)/2 and ry>(yr-yl)/2
    return A.bounded((xr-xl)*(yr-yl)),((xr-xl)/(2*rx),(yr-yl)/(2*ry))


def choose(area,moduli,ratios,orders,allowance):
    area=A.bounded(Q(area));allowance=A.bounded(Q(allowance))
    assert area>0 and allowance>0 and len(moduli)==len(ratios)==2
    moduli=tuple(A.bounded(Q(q)) for q in moduli)
    ratios=tuple(A.bounded(Q(q)) for q in ratios)
    assert all(q>=0 for q in moduli) and all(0<q<1 for q in ratios)
    orders=tuple(orders);assert orders and orders==tuple(sorted(set(orders)))
    assert all(type(n) is int and 2<=n<=64 for n in orders)
    target=A.bounded(3*allowance/4)
    for n in orders:
        error,parts=RES.error(area,moduli,ratios,(n,n))
        if error<=target:
            return {'order':n,'tensor_nodes':n*n,'remainder_bound':error,
                'axis_remainders':parts,'original_allowance':allowance,'remainder_budget':target,
                'reserved_central_radius':A.bounded(allowance-error)}
    return None


def audit_choice(choice,area,moduli,ratios,orders,allowance):
    assert choice==choose(area,moduli,ratios,orders,allowance),'Changed uniform-order choice'
    assert choice is not None
    return choice


def enclosure(data,N,n,box,radii,allowance):
    area,ratios=bounds(box,radii)
    # This also verifies every original majorant algebra witness. It does not
    # require the historical broad final interval to fit the new cell budget.
    A.cell(data,N,n,area,*ratios)
    assert data['method']=='tensor_Gauss_Cauchy' and data['central']['finite']
    moduli=[A.interval(data[k])[1] for k in ('majorant_x','majorant_y')]
    error,parts=RES.error(area,moduli,ratios,(n,n))
    coarse=2*area*sum(m*q**(2*n)/(1-q) for m,q in zip(moduli,ratios))
    assert error<=coarse,'The proven bound must improve the historical Taylor bound'
    centre=A.interval(data['central']['real']);imag=A.interval(data['central']['imag'])
    assert imag[0]<=0<=imag[1]
    interval=(A.bounded(centre[0]-error),A.bounded(centre[1]+error))
    radius=A.bounded((interval[1]-interval[0])/2);allowance=A.bounded(Q(allowance))
    assert allowance>0 and radius<=allowance,'Sharper certificate still exceeds original allocation'
    return {'interval':interval,'central':centre,'error':error,'axis_remainders':parts,
        'original_allowance':allowance,'radius':radius,'area':area,'ratios':ratios,
        'order':n,'tensor_nodes':n*n,'independent_geometry_and_central_replay_required':True,
        'whole_graph_certificate':False}
