"""Exact homogeneous real-box bound, independently rebuilt from coefficients."""
from pathlib import Path
from fractions import Fraction as Q
import theta_two_radial_polynomial_audit as B
import exact_exp2_enclosures as E
import exact_radial_node_bounds_v2 as L
import theta_two_radial_majorant_audit_v4 as A


class Homogeneous(B.Polynomial):
    def __init__(self,path):
        with Path(path).open() as stream:
            self.order,count=map(int,stream.readline().split());assert 1<=self.order<=22 and count<=100000
            degree=12*self.order-2;terms={}
            for line in stream:
                i,j,n,d=map(int,line.split());assert i>=0 and j>=0 and i+j<=degree and d>0
                key=degree-i-j,j;assert key not in terms;terms[key]=B.scaled(abs(Q(n,d)))
        assert len(terms)==count;self.count=count;self.rows=[]
        for i in range(max(i for i,j in terms)+1):
            js=[j for ii,j in terms if ii==i]
            self.rows.append([terms.get((i,j),(0,0)) for j in range(max(js)+1)] if js else [])

    def verify_homogeneous(self,proof):
        return super().verify({'r_mod':proof['s_hi'],'w_mod':proof['v_hi'],
                               'numerator_absolute':proof['numerator_normalized']})


def verify(data,request,polynomial):
    words=request.split();assert words[0]=='normalized-range' and len(words)==11
    bits,chart=map(int,words[1:3]);assert 256<=bits<=512 and chart in (0,1)
    xl,xr,yl,yr=[Q(int(words[j]),int(words[j+1])) for j in range(3,11,2)]
    assert -128<=xl<xr<=128 and 0<=yl<yr<=(128 if chart else 1)
    assert data['finite'] and data['method']=='normalized_radial_range' and data['chart']==chart
    N=polynomial.order;assert data['maximum_order']==N
    area=(xr-xl)*(yr-yl);assert A.interval(data['area'])[0]<=area<=A.interval(data['area'])[1]
    p={}
    for k,v in data['homogeneous_majorant'].items():
        lo,hi=A.interval(v);assert lo==hi and lo>=0;p[k]=lo
    rl,rh,wl,wh=(p[k] for k in ('r_lo','r_hi','w_lo','w_hi'))
    assert 0<rl<=E.exp2(xl)[0] and rh>=E.exp2(xr)[1]
    assert wl<=(E.exp2(yl)[0] if chart else yl) and wh>=(E.exp2(yr)[1] if chart else yr)
    sl,sh,vl,vh=(p[k] for k in ('s_lo','s_hi','v_lo','v_hi'))
    assert 0<sl<=1/rh and sh>=1/rl and vl<=wl/rh and vh>=wh/rl
    assert 0<p['T_lo']<=(1+2*sl)*(1+2*vl)/4 and p['T_hi']>=(1+2*sh)*(1+2*vh)/4
    assert 0<p['U_lo']<=(1+2*vl+2*sl)/4 and p['U_hi']>=(1+2*vh+2*sh)/4
    assert 0<p['Lt']<=sl*sl+p['T_lo'] and 0<p['Lu']<=sl*sl+p['U_lo']
    assert 0<p['D_lo']<=1+vl+sl and 0<p['AD_lo']<=(1+2*sl)*(1+2*vl+2*sl)/2
    assert 2<=p['t_lower']<=rl*rl*p['T_lo'] and 2<=p['u_lower']<=rl*rl*p['U_lo']
    assert p['log_t']>=L.log_bounds(rh*rh*p['T_hi'])[1]
    assert p['log_u']>=L.log_bounds(rh*rh*p['U_hi'])[1]
    assert p['ln2_upper']>=E.log2()[1]
    radial=vh*p['ln2_upper']**2/rl**2 if chart else p['ln2_upper']/rl**3
    assert p['radial_scale']>=radial
    coefficient_replay=polynomial.verify_homogeneous(data['homogeneous_majorant'])
    need=p['radial_scale']*p['T_hi']**3*p['U_hi']**3*p['log_t']*p['log_u']*p['numerator_normalized']
    need/=144*p['Lt']**(2*N+3)*p['Lu']**(2*N+3)*p['D_lo']**(4*N-2)*p['AD_lo']**2
    assert p['bound']>=need
    lo,hi=A.interval(data['real']);assert lo<=-p['bound']*area and hi>=p['bound']*area
    return {'upper_density':p['bound'],'interval':(lo,hi),'area':area,
            'coefficient_replay':coefficient_replay,'scalar_and_transcendental_premises_exact':True}
