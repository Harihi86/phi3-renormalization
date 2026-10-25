"""Exact rational map followed by the independent original routed reference."""
import exact_histogram_transport_v1 as A
import g371_histogram_models_v1 as P
import g371_soft_reference_v1 as H
M=H.M
def evaluate(point,binding,digits=640,mode='hist-pair',quadrature=False):
 assert mode in ('hist-pair','hist-scale')
 base='soft-pair' if mode=='hist-pair' else 'soft-scale'
 mapped,jac=A.transform(P.MODELS[base.replace('-','_')],point)
 ctx,value=H.evaluate(mapped,binding,digits,base,quadrature)
 return ctx,value*ctx.mpf(jac.numerator)/jac.denominator
