"""Randomized QMC of the complete graph-compiled BPH recursive density.

Importance proposals depend on graph/program structure, never on graph IDs.
The defensive recursive mixture weights every point by its complete PDF.
"""
from math import exp, isfinite, lgamma, sqrt
from time import perf_counter
from .automatic_bph import compile_bph


def integrate_compiled(compiled,invariants=None,m2=1,*,power=14,replicates=8,seed=20261019,beta=.35,chunk_size=2048,sampler='auto'):
    try:
        import numpy as np
        from scipy.special import gammaincinv
        from scipy.stats import qmc
    except ImportError as exc:
        raise RuntimeError('Automatic integration requires optional NumPy and SciPy; expression compilation itself uses only the standard library') from exc
    if not hasattr(compiled,'batch_density'):compiled=compile_bph(compiled)
    if type(power) is not int or not 6<=power<=22:raise ValueError('QMC power must be an integer from6 to22')
    if type(replicates) is not int or not 2<=replicates<=32:raise ValueError('Use2--32 independent randomized replicates')
    if not isfinite(float(beta)) or not .15<=beta<=1:raise ValueError('Dirichlet beta must lie between.15 and1')
    if type(chunk_size) is not int or chunk_size<1:raise ValueError('Positive integer chunk size is required')
    if sampler not in ('auto','hepp','dirichlet','recursive','mixture'):raise ValueError('Sampler must be auto, hepp, dirichlet, recursive or mixture')
    if sampler=='auto':sampler='hepp' if not compiled.parts else 'mixture'
    start=perf_counter();g=compiled.graph;size=len(g.edges);normalization=exp(size*lgamma(beta)-lgamma(size*beta))
    estimates=[];half_estimates=[];maximum=0.;negative=0;precision={}
    for repeat in range(replicates):
        dimension=(1+max(2*(size-1),size) if sampler=='mixture' else
                   2*(size-1) if sampler in ('hepp','recursive') else size)
        points=qmc.Sobol(dimension,scramble=True,seed=seed+104729*repeat).random_base2(power)
        full_sum=0.;half_sum=0.;offset=0
        for begin in range(0,len(points),chunk_size):
            uniforms=points[begin:begin+chunk_size]
            uniforms=np.clip(uniforms,np.nextafter(0.,1.),np.nextafter(1.,0.))
            if sampler in ('recursive','mixture'):
                from .recursive_sampling import recursive_hepp_samples,recursive_mixture_samples
                sample=recursive_mixture_samples if sampler=='mixture' else recursive_hepp_samples
                alpha,weights,sampling_metadata=sample(compiled,uniforms)
            elif sampler=='hepp':
                from .automatic_sampling import hepp_samples
                alpha,weights,sampling_metadata=hepp_samples(g,uniforms)
            else:
                gamma=gammaincinv(beta,uniforms)
                alpha=gamma/gamma.sum(axis=1,keepdims=True)
                weights=normalization*np.exp((1-beta)*np.log(alpha).sum(axis=1))
                sampling_metadata={'dirichlet_beta':beta}
            density=compiled.batch_density(alpha,invariants,m2,diagnostics=precision)
            values=density*weights
            if not np.all(np.isfinite(values)):raise ArithmeticError('Nonfinite QMC weighted density; use another distribution or a resolved chart')
            full_sum+=float(np.sum(values,dtype=np.longdouble))
            take=max(0,min(len(values),2**(power-1)-begin))
            half_sum+=float(np.sum(values[:take],dtype=np.longdouble))
            maximum=max(maximum,float(np.max(np.abs(values))));negative+=int(np.sum(values<0));offset+=len(values)
        estimates.append(full_sum/offset);half_estimates.append(half_sum/2**(power-1))
    mean=float(np.mean(estimates));stderr=float(np.std(estimates,ddof=1)/sqrt(replicates))
    return {'value':mean,'standard_error':stderr,'replicate_values':estimates,
            'half_power_replicate_values':half_estimates,
            'half_power_value':float(np.mean(half_estimates)),
            'power_refinement_difference':abs(mean-float(np.mean(half_estimates))),
            'points_per_replicate':2**power,'replicates':replicates,'total_samples':replicates*2**power,
            'seed':seed,'sampler':sampler,'sampling_metadata':sampling_metadata,'maximum_absolute_weighted_sample':maximum,
            'negative_weighted_samples':negative,'seconds':perf_counter()-start,
            'precision_control':precision,
            'method':'Complete compiled recursive Gaussian BPH sum; independently scrambled Sobol with graph-derived importance sampling',
            'topology_specific_kernel_used':False,
            'error_scope':'Between-scramble standard error and nested-power difference are numerical diagnostics, not rigorous quadrature enclosures'}
