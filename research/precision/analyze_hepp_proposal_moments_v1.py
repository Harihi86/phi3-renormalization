"""Exact proposal-only moments. No physics evaluation or acceptance decision."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
from collections import Counter
from datetime import datetime
import hashlib,json,math

HERE=Path(__file__).resolve().parent
TABLES=HERE/'evidence/named-compact-20261025-101856/mixed-forest-hepp-proof-01/tables.json'
TABLE_SHA='83e7e3bbc4a3ca36a84fdb8eb52f87d4b51cb8fdd2342a9ba28c85120a67b7d1'


def chain_sum(factors,n=8):
    """Sum products of proper-subset factors over every edge order."""
    full=(1<<n)-1;values={0:Q(1)}
    for size in range(1,n):
        for mask in range(1,full):
            if mask.bit_count()==size:
                values[mask]=factors[mask]*sum((values[mask^(1<<e)] for e in range(n) if mask&(1<<e)),Q())
    return sum((values[full^(1<<e)] for e in range(n)),Q())


def mobius_second(beta,c):
    """Integral of r^(2 beta-2)/h(r), h=c/[c-(c-1)r]^2."""
    assert beta>=1 and c>0
    return (c*c/(2*beta-1)-2*c*(c-1)/(2*beta)+(c-1)**2/(2*beta+1))/c


def encode(x):
    if isinstance(x,Q):return {'numerator':str(x.numerator),'denominator':str(x.denominator)}
    if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def main():
    assert hashlib.sha256(TABLES.read_bytes()).hexdigest()==TABLE_SHA
    # Independent explicit enumeration for the subset recursion.
    factors={s:Q(s+3,2*s+5) for s in range(1,15)}
    direct=Q()
    for order in permutations(range(4)):
        remaining=15;term=Q(1)
        for edge in order[:-1]:remaining^=1<<edge;term*=factors[remaining]
        direct+=term
    assert direct==chain_sum(factors,4)
    assert chain_sum({s:Q(1,s.bit_count()) for s in range(1,255)})==8
    # The rational inverse map has exact endpoints and positive derivative.
    # h(r) follows by solving u=r/[c-(c-1)r].
    identities=0
    for beta in map(Q,range(1,9)):
        c=beta
        for u in (Q(),Q(1,7),Q(1,2),Q(5,6),Q(1)):
            r=c*u/(1+(c-1)*u)
            inverse=r/(c-(c-1)*r)
            derivative=c/(1+(c-1)*u)**2
            density=c/(c-(c-1)*r)**2
            assert inverse==u and derivative*density==1 and derivative>0
            identities+=1
        assert mobius_second(beta,c)>=1/beta**2
        assert mobius_second(beta,c)<=mobius_second(beta,Q(1))
    tables=json.loads(TABLES.read_text());result={}
    for gid,variants in tables.items():
        result[gid]={}
        for name,table in variants.items():
            powers={int(s):p for s,p in table['powers'].items()}
            gamma={int(s):int(v) for s,v in table['gamma'].items()}
            # This clipped monomial proxy is not the signed physical density.
            beta={s:max(Q(1,2),Q(gamma[s])) for s in powers}
            normal=chain_sum({s:Q(p) for s,p in powers.items()})
            assert normal==table['normalization']
            constant_second=normal*chain_sum({s:Q(p,2*p*s.bit_count()-1) for s,p in powers.items()})
            proxy_mean=chain_sum({s:1/b for s,b in beta.items()})
            current_second=normal*chain_sum({s:Q(p)/(2*p*beta[s]-1) for s,p in powers.items()})
            # Proposed order mass prod(1/beta)/Z and a rational ratio map:
            # r=u^2 for beta=1/2, else r=beta*u/[1+(beta-1)u].
            factors={s:(1/b if b==Q(1,2) else b*mobius_second(b,b)) for s,b in beta.items()}
            candidate_second=proxy_mean*chain_sum(factors)
            assert constant_second>=64 and current_second>=proxy_mean**2 and candidate_second>=proxy_mean**2
            result[gid][name]={
                'gamma_histogram':dict(sorted(Counter(gamma.values()).items())),
                'constant_integrand_mean':Q(8),
                'constant_single_selector_second_over_mean_squared':constant_second/64,
                'proxy_mean':proxy_mean,
                'proxy_current_single_selector_second_over_mean_squared':current_second/proxy_mean**2,
                'proxy_mobius_single_selector_second_over_mean_squared':candidate_second/proxy_mean**2,
                'proxy_variance_candidate_over_current':(candidate_second-proxy_mean**2)/(current_second-proxy_mean**2),
                'floats_for_display':{'constant_relative_second_moment':float(constant_second/64),
                    'proxy_current_relative_second_moment':float(current_second/proxy_mean**2),
                    'proxy_candidate_relative_second_moment':float(candidate_second/proxy_mean**2)}}
    output={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'tables_sha256':TABLE_SHA,'exact_map_identities':identities,'exact_recursion_vs_explicit_orders':24,
        'new_physics_calls':0,'physical_acceptances_added':0,
        'scope':'Exact single-selector moments of a constant and a clipped sector monomial proxy; neither is the actual signed forest density. Eight-rotation averaging has no larger second moment by Jensen, but its actual covariance is not calculated here. No precision gate or production policy is changed.',
        'proxy_definition':'On an ordered sector, f_proxy=product_j r_j^(beta(S_j)-|S_j|), beta=max(gamma,1/2).',
        'result':result}
    path=HERE/'evidence'/('hepp-proposal-moments-'+datetime.now().strftime('%Y%m%d-%H%M%S')+'.json')
    with path.open('x',encoding='utf-8') as handle:json.dump(encode(output),handle,indent=2)
    print(json.dumps({'receipt':str(path),'exact_map_identities':identities,'new_physics_calls':0,
        'results':{g:{v:r['floats_for_display'] for v,r in rows.items()} for g,rows in result.items()}}))

if __name__=='__main__':main()
