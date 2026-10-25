"""Overall T0 common-scale primitives for shared or distinct bubble edges."""
import double_weighted_e2_scalar_v1 as D
import shared_weighted_e2_scalar_v1 as H

def vertex_J(ctx,mass,edge,z):
    m,d,z=map(ctx.mpf,(mass,edge,z));assert m>d>0 and z>=0
    return -D.F.M.bubble_primitive(ctx,(m+z-d)/d)+D.F.M.bubble_primitive(ctx,(m-d)/d)

def vertex_J_mass(ctx,mass,edge,z):
    m,d,z=map(ctx.mpf,(mass,edge,z));a=(m-d)/d;b=(m+z-d)/d
    return -(D.B.bubble(ctx,b)/(1+b)**2-D.B.bubble(ctx,a)/(1+a)**2)/d

def distinct_kernel(ctx,mass,first,second,z,tau):
    m,d,e,z,t=map(ctx.mpf,(mass,first,second,z,tau));assert m>d+e and min(d,e)>0
    c,h,k=D.coefficients(ctx,t)
    return c*vertex_J(ctx,m,d,z)-h*e*vertex_J_mass(ctx,m,d,z)+k/t*vertex_J(ctx,m+e*(1/t-1),d,z)

def nine_term_reference(ctx,mass,first,second,z,tau):
    m,d,e,z,t=map(ctx.mpf,(mass,first,second,z,tau));ct,dt,et=D.coefficients(ctx,t);shift=e*(1/t-1)
    def R(v):return -ctx.log1p(z/v)
    def Rm(v):return z/(v*(v+z))
    def Rmm(v):return -z*(2*v+z)/(v*v*(v+z)**2)
    def integrand(u):
        s=(1-u*u)/4
        if not s:
            def limit(v):return -2*R(v)-d*Rm(v)-z/d
            def limit_m(v):return -2*Rm(v)-d*Rmm(v)
            total=ct*limit(m)-dt*e*limit_m(m)+et/t*limit(m+shift)
        else:
            cs,ds,es=D.coefficients(ctx,s);shift_s=d*(1/s-1)
            total=cs*ct*R(m)-cs*dt*e*Rm(m)+cs*et/t*R(m+shift)
            total+=-ds*d*ct*Rm(m)+ds*d*dt*e*Rmm(m)-ds*d*et/t*Rm(m+shift)
            total+=es/s*ct*R(m+shift_s)-es/s*dt*e*Rm(m+shift_s)+es/s*et/t*R(m+shift_s+shift)
        return u**4/12*total
    return ctx.quad(integrand,[0,ctx.mpf('.5'),1])

def shared_density(ctx,mass,edge,z,v,bubble=None):
    m,d,z,v=map(ctx.mpf,(mass,edge,z,v));assert m>d>0 and z>=0 and 0<=v<=1
    b=D.B.bubble if bubble is None else bubble;a=m+v*z
    return -d*d*z*b(ctx,(a-d)/d)**2/a**3

def shared_spectral_primitive(ctx,mass,edge,z,s,t):
    m,d,z,s,t=map(ctx.mpf,(mass,edge,z,s,t));c1,c2,c3,cs,ct=H.partial_coefficients(ctx,s,t)
    def r(a):return -ctx.log1p(z/a)
    rm=z/(m*(m+z));rmm=-z*(2*m+z)/(m*m*(m+z)**2)
    return c1*r(m)-c2*d*rm+c3*d*d*rmm/2+cs/s*r(m+d*(1/s-1))+ct/t*r(m+d*(1/t-1))

def shared_spectral_reference(ctx,mass,edge,z,s,t):
    m,d,z,s,t=map(ctx.mpf,(mass,edge,z,s,t))
    def f(v):
        x=(m+v*z-d)/d
        return -x**4/((1+x)**3*(1+s*x)*(1+t*x))/d
    return z*ctx.quad(f,[0,ctx.mpf('.5'),1])
