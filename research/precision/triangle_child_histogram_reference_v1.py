"""Exact rational map, followed by independent original-route scalar AD reference."""
import exact_histogram_transport_v1 as A
import triangle_child_histogram_models_v1 as P
import triangle_two_host_reference_v1 as H
M=H.M
AD=H.AD
F=H.F
N=H.N
def evaluate(point,binding,digits=640,mode='hist-pair',quadrature=False):
 assert mode in ('hist-pair','hist-scale')
 base='triangle-child-pair' if mode=='hist-pair' else 'triangle-child-scale'
 mapped,jac=A.transform(P.MODELS[binding['id']][base.replace('-','_')],point)
 ctx,value=H.evaluate(mapped,binding,digits,base,quadrature)
 return ctx,value*ctx.mpf(jac.numerator)/jac.denominator
