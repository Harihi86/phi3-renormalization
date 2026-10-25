"""Exact two-sector radius-difference transport of the complete G20 q-angle integral."""
import double_weighted_loop_angle_v2 as M


def evaluate(ctx,point,chart='linear',terms=1024,chebyshev=False):
    x,y,u=map(ctx.mpf,point);assert all(0<v<1 for v in (x,y,u))
    power={'linear':1,'cubic':3}[chart]
    convert=lambda v:(v/(1-v)**power,(1+(power-1)*v)/(1-v)**(power+1))
    radius,dr=convert(x);difference,dd=convert(y)
    t=radius*radius;s=(radius+difference)**2;jac=4*radius*(radius+difference)*dr*dd
    value=M.angular_orbit(ctx,s,t,u,8,terms,chebyshev)+M.angular_orbit(ctx,t,s,u,8,terms,chebyshev)
    return jac*s*s*t*t*M.G.B.bubble(ctx,s)*M.G.B.bubble(ctx,t)/(4*(1+s)**2*(1+t)**2)*value
