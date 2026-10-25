"""Exact interval, endpoint and API checks for both signed-density components."""
from fractions import Fraction as Q
from fourier_spectral_audit import interval,exp_negative64

def verify(data):
    assert data['success'] and 192<=data['bits']<=1024 and 16<=data['evaluation_cap']<=20000
    rl,rh=interval(data['r']);assert Q(1,262144)<=rl<=rh<=12
    components=data['components'];assert len(components)==2
    el,eh=exp_negative64();rows=[]
    for row,name,constant in zip(components,('plus','minus'),(2000,200)):
        assert row['name']==name and row['passed'] and row['calc_exit']==0 and row['calls']>=row['analytic_calls']>0
        ql,qh=interval(row['quadrature']);il,ih=interval(row['quadrature_imag']);assert il<=0<=ih
        tl,th=interval(row['tail']);assert tl<=constant*el and th>=constant*eh
        cl,ch=interval(row['C_scaled']);assert 0<cl<=max(Q(),ql) and ch>=qh+th
        assert (ch-cl)/2<(1+abs((ch+cl)/2))/10**12
        rows.append({'name':name,'C_scaled':(cl,ch),'quadrature':(ql,qh),'tail_upper':th})
    return {'r_interval':(rl,rh),'components':rows,'FLINT_is_interval_bound_authority':True,'positive_components_not_a_positive_signed_density':True}
