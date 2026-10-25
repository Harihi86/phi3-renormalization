"""Exact stored-mesh audit, independent of native interval arithmetic."""
from fractions import Fraction as Q
def interval(v):
    a,b,e=map(int,v);assert a<=b;return Q(a)*Q(2)**e,Q(b)*Q(2)**e
def exp_negative64():
    term=Q(1);total=term
    for n in range(1,385):term*=Q(64,n);total+=term
    nxt=term*Q(64,385);upper=total+nxt/(1-Q(64,386))
    return 1/upper,1/total
def verify(data):
    assert data['success'] and data['invalid_cells']==0 and data['bits']>=128
    rlo,rhi=interval(data['r']);assert Q(1,262144)<=rlo<=rhi<=12
    leaves={(0,0)}
    for i,l in data['splits']:
        assert (i,l) in leaves and l<30;leaves.remove((i,l));leaves.update(((2*i,l+1),(2*i+1,l+1)))
    assert len(leaves)==len(data['cells'])<=data['cell_cap'] and {tuple(c['geometry']) for c in data['cells']}==leaves
    sl=sh=eh=Q()
    for c in data['cells']:
        assert c['valid'];i,l=c['geometry'];a,b=Q(32*i,2**l),Q(32*(i+1),2**l)
        assert interval(c['left'])==(a,a) and interval(c['right'])==(b,b)
        dl,dh=interval(c['sixth']);need=(b-a)**7*max(abs(dl),abs(dh))/2800
        el,ehi=interval(c['error']);assert 0<=el<=ehi and ehi>=need
        cl,ch=interval(c['center']);sl+=cl;sh+=ch;eh+=ehi
    el,eu=exp_negative64();tl,tu=interval(data['spectral_tail']);assert tl<=Q(8,27)*396*el and tu>=Q(8,27)*396*eu
    lo,hi=interval(data['C_scaled']);assert lo<=max(Q(),sl-eh) and hi>=sh+eh+tu
    assert (hi-lo)/2<(1+abs((hi+lo)/2))/10**12
    return {'cells':len(leaves),'quadrature_center_interval':(sl,sh),'quadrature_error_upper':eh,'spectral_tail_upper':tu,'C_scaled':(lo,hi)}
