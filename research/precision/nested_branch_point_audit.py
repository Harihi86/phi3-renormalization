"""Independent rational reconstruction of positive nested spectral enclosures."""
from fractions import Fraction as Q
from functools import lru_cache
from fourier_spectral_audit import interval,exp_negative64

@lru_cache(maxsize=1)
def exp_negative288():
    term=total=Q(1)
    for n in range(1,1153):term*=Q(288,n);total+=term
    upper=total+term*Q(288,1153)/(1-Q(288,1154))
    return 1/upper,1/total

def verify(data):
    assert data['success'] and data['calc_exit']==0
    assert 256<=data['bits']<=1024 and 16<=data['evaluation_cap']<=20000
    assert not data['exhausted'] and data['inner_failures']==0 and data['calls']>=data['analytic_calls']>0
    rl,rh=interval(data['r']);assert Q(1,262144)<=rl<=rh<=12
    if data['component']=='continuum':
        assert 0<data['inner_calls']<=500000 and data['wide_calls']>0 and data['narrow_calls']>0
        lo,hi=exp_negative288();factor=2**96
    else:
        assert data['component']=='mass' and data['inner_calls']==data['wide_calls']==data['narrow_calls']==0
        lo,hi=exp_negative64();factor=50
    ql,qh=interval(data['quadrature']);il,ih=interval(data['quadrature_imag']);assert il<=0<=ih
    tl,th=interval(data['tail']);assert tl<=factor*lo and th>=factor*hi
    cl,ch=interval(data['C_scaled']);assert 0<cl<=max(Q(),ql) and ch>=qh+th
    assert (ch-cl)/2<(1+abs((ch+cl)/2))/10**12
    return {'r_interval':(rl,rh),'component':data['component'],'C_scaled':(cl,ch),'tail_upper':th,'FLINT_is_interval_bound_authority':True}
