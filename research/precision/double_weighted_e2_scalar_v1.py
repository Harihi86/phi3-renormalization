"""One remaining spectral parameter for two distinct renormalized bubble edges."""
import weighted_bubble_e2_closed_v1 as F
import weighted_bubble_e2_integral_reference_v1 as B

def mass_derivative(ctx,mass,edge,z):
    m,d,z=map(ctx.mpf,(mass,edge,z));assert m>d>0 and z>=0
    t0=(m-d)/d;t1=t0+z/d
    return F.M.bubble_primitive(ctx,t1)-F.M.bubble_primitive(ctx,t0)-(z/d)*B.bubble(ctx,t0)/(1+t0)**2

def coefficients(ctx,sigma):
    s=ctx.mpf(sigma);assert 0<s<=ctx.mpf(1)/4
    return (s-2)/(1-s)**2,1/(1-s),1/(1-s)**2

def spectral_kernel(ctx,mass,first,second,z,tau):
    m,d,e,z,t=map(ctx.mpf,(mass,first,second,z,tau));assert m>d+e and min(d,e)>0
    c,h,k=coefficients(ctx,t)
    return c*F.closed_J(ctx,m,d,z)-h*e*mass_derivative(ctx,m,d,z)+k/t*F.closed_J(ctx,m+e*(1/t-1),d,z)

def nine_term_reference(ctx,mass,first,second,z,tau):
    """Integrate the first spectral variable using elementary nine-term expansion only."""
    m,d,e,z,t=map(ctx.mpf,(mass,first,second,z,tau));ct,dt,et=coefficients(ctx,t);shift=e*(1/t-1)
    def R(v):return (v+z)*ctx.log1p(z/v)-z
    def Rm(v):return ctx.log1p(z/v)-z/v
    def Rmm(v):return z*z/(v*v*(v+z))
    def integrand(u):
        s=(1-u*u)/4
        if not s:
            def limit(v):return -2*R(v)-d*Rm(v)+z*z/(2*d)
            def limit_m(v):return -2*Rm(v)-d*Rmm(v)
            total=ct*limit(m)-dt*e*limit_m(m)+et/t*limit(m+shift)
        else:
            cs,ds,es=coefficients(ctx,s);shift_s=d*(1/s-1)
            total=cs*ct*R(m)-cs*dt*e*Rm(m)+cs*et/t*R(m+shift)
            total+=-ds*d*ct*Rm(m)+ds*d*dt*e*Rmm(m)-ds*d*et/t*Rm(m+shift)
            total+=es/s*ct*R(m+shift_s)-es/s*dt*e*Rm(m+shift_s)+es/s*et/t*R(m+shift_s+shift)
        return u**4/12*total
    return ctx.quad(integrand,[0,ctx.mpf('.5'),1])
