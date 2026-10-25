"""Exact rational map, followed by independent original-route scalar AD reference."""
import exact_histogram_transport_v1 as A
import bubble_child_histogram_models_v1 as P
import bubble_two_host_reference_v2 as H
M=H.M
AD=H.AD
KA=H.KA
N=H.N
def evaluate(point,binding,digits=640,mode='hist-pair',quadrature=False):
 assert mode in ('hist-pair','hist-scale')
 base='child-pair' if mode=='hist-pair' else 'child-scale'
 mapped,jac=A.transform(P.MODELS[binding['id']][base.replace('-','_')],point)
 ctx,value=H.evaluate(mapped,binding,digits,base,quadrature)
 return ctx,value*ctx.mpf(jac.numerator)/jac.denominator
