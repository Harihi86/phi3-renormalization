"""Exact Fraction transport followed by the independent original Gaussian Schur reference."""
import mixed_single_weighted_e2_reference_v1 as B
import exact_histogram_transport_v1 as A
import mixed_schwinger_histogram_models_v1 as P
J,X=B.J,B.X
def evaluate(points,binding,digits=640,chart='hist-linear',quadrature=False):
    assert chart in ('hist-linear','hist-square')
    base=chart.replace('hist','weighted');model=P.MODELS[binding['id']][base.replace('-','_')]
    mapped,jac=A.transform(model,points)
    ctx,result=B.evaluate(mapped,binding,digits,base,quadrature)
    return ctx,result*ctx.mpf(jac.numerator)/jac.denominator
