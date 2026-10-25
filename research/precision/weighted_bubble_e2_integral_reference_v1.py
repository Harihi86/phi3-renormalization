"""Independent full remainder integral using an elementary bubble expression."""
def bubble(ctx,t):
    t=ctx.mpf(t)
    if not t:return ctx.mpf(0)
    return (t+4)*ctx.sqrt((t+4)/t)*ctx.asinh(ctx.sqrt(t)/2)/3-4*t/9-ctx.mpf(4)/3

def reference_J(ctx,mass,edge,z):
    mass,edge,z=map(ctx.mpf,(mass,edge,z))
    assert mass>edge>0 and z>=0
    if not z:return ctx.mpf(0)
    # Add a scale-aware division point when z is large. This resolves the
    # narrow transition without dropping any part of the integration domain.
    cut=min(ctx.mpf('0.5'),mass/z)
    return edge*z*z*ctx.quad(lambda u:(1-u)*bubble(ctx,(mass+u*z-edge)/edge)/(mass+u*z)**2,[0,cut,ctx.mpf('0.5'),1])
