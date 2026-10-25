"""Exact unchanged four-replicate gates for complete transformed integrals."""
from fractions import Fraction as Q
from decimal import Decimal,localcontext


def functionals(bounds):
    count=len(bounds)
    assert count>=2 and count&(count-1)==0
    assert all(type(a) is Q and type(b) is Q and a<=b for a,b in bounds)
    def mean(rows):
        lo=sum((r[0] for r in rows),Q())/len(rows)
        hi=sum((r[1] for r in rows),Q())/len(rows)
        return {'midpoint':(lo+hi)/2,'radius':(hi-lo)/2}
    full,half=mean(bounds),mean(bounds[:count//2])
    gates={k:v['radius']<=Q.from_float(2e-12)*abs(v['midpoint']) for k,v in (('full',full),('half',half))}
    return {'full':full,'half':half,'nodes':count,'arithmetic_components':gates,'arithmetic_passed':all(gates.values())}


def statistics(groups):
    assert set(groups)=={'atan','asin'} and all(len(rows)==4 for rows in groups.values())
    output={}
    for name,rows in groups.items():
        assert all(row['arithmetic_passed'] for row in rows)
        values=[r['full']['midpoint'] for r in rows];halves=[r['half']['midpoint'] for r in rows]
        mean=sum(values,Q())/4;half=sum(halves,Q())/4
        variance=sum(((v-mean)**2 for v in values),Q())/12
        delta=abs(mean-half);target=Q(1,10**8)+Q(3,100)*abs(mean)
        output[name]={'replicate_means':values,'half_replicate_means':halves,'mean':mean,'half_mean':half,
            'standard_error_squared':variance,'half_difference':delta,'own_target':target,
            'own_passed':variance<=target*target and delta<=target}
    g,p=output['atan'],output['asin'];centre=(g['mean']+p['mean'])/2
    target=Q(1,10**8)+Q(3,100)*abs(centre);gap=abs(g['mean']-p['mean'])
    variance=g['standard_error_squared']+p['standard_error_squared']
    halves=2*(g['half_difference']+p['half_difference'])
    conditions={'root_own':g['own_passed'],'center_own':p['own_passed'],
        'pooled_centre_gap':gap<=target,'four_SE_plus_two_halves':target>=halves and 16*variance<=(target-halves)**2}
    with localcontext() as ctx:
        ctx.prec=70
        dec=lambda q:Decimal(q.numerator)/Decimal(q.denominator)
        diagnostic={'pooled_mean':str(dec(centre)),'pooled_target':str(dec(target)),
            'centre_gap_ratio':str(dec(gap)/dec(target)),
            'strong_uncertainty_ratio':str((4*dec(variance).sqrt()+dec(halves))/dec(target))}
    return {'kernels':output,'pooled_centre':centre,'pooled_target':target,'absolute_centre_gap':gap,
        'combined_SE_squared':variance,'two_half_differences':halves,'conditions':conditions,
        'completion_gate_passed':all(conditions.values()),'diagnostics':diagnostic,
        'rigorous_quadrature_error_bound':False}
