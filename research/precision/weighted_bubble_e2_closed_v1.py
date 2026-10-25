"""An exact closed primitive for B(t)/(1+t)^2 in an overall T2 subtraction."""
import nested_analytic_marginals_v1 as M
import overlapping_bubble_self_energy_factor_v1 as O

def weighted_moment(ctx,t):
    t=ctx.mpf(t)
    if not t:return ctx.mpf(0)
    h=ctx.asinh(ctx.sqrt(t)/2)
    return h*ctx.sqrt(t*(t+4))*(t-2)/(3*(t+1))+4*h*h/3-11*t/18+ctx.log1p(t)/18+8*t/(9*(t+1))

def closed_J(ctx,mass,edge,z):
    mass,edge,z=map(ctx.mpf,(mass,edge,z));assert mass>edge>0 and z>=0
    if not z:return ctx.mpf(0)
    t0=(mass-edge)/edge;t1=t0+z/edge
    return edge*(t1*(M.bubble_primitive(ctx,t1)-M.bubble_primitive(ctx,t0))-(weighted_moment(ctx,t1)-weighted_moment(ctx,t0)))

def reference_J(ctx,mass,edge,z):return O.primitive_J(ctx,mass,edge,z)
