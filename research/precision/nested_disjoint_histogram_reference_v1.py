"""Exact rational map, followed by independent original-route scalar AD reference."""
import exact_histogram_transport_v1 as A
import nested_disjoint_histogram_models_v1 as P
import nested_disjoint_reference_v1 as H
M=H.M
AD=H.AD
def evaluate(point,binding,digits=640,mode='hist-pair',quadrature=False):
 assert mode in ('hist-pair','hist-scale')
 base='nested-pair' if mode=='hist-pair' else 'nested-scale'
 mapped,jac=A.transform(P.MODELS[binding['id']][base.replace('-','_')],point)
 ctx,value=H.evaluate(mapped,binding,digits,base,quadrature)
 return ctx,value*ctx.mpf(jac.numerator)/jac.denominator
