"""Exact one-dimensional conditional primitive, shared only as mathematics."""


def mean_h(mp,t,M,A,B,C):
    """Integral_0^1 h(t*(-A*x*x+B*x+C)/M) dx; all inputs positive real."""
    t,M,A,B,C=map(mp.mpf,(t,M,A,B,C))
    h=lambda z:(1+z)*mp.log1p(z)-z
    if not A:
        assert not B
        return h(t*C/M)
    centre=B/(2*A);D=t*A/M;S=1+t*(C+B*B/(4*A))/M
    def primitive(v):
        f=S-D*v*v
        assert f>0
        return (S*v-D*v**3/3)*mp.log(f)+5*D*v**3/9+(1-7*S/3)*v+4*S**mp.mpf('1.5')/(3*mp.sqrt(D))*mp.atanh(v*mp.sqrt(D/S))
    return primitive(1-centre)-primitive(-centre)


def core(mp,t,a,b,y):
    t,a,b,y=map(mp.mpf,(t,a,b,y));U=a*b+a+b;M=1+a+b
    A=a*a*(b+1)/U;B=a-2*a*b*y/U;C=b*y*(1-y)+a*b*y*y/U
    return a*b*M*(mean_h(mp,t,M,A,B,C)/U**3
          -mean_h(mp,t,M,0,0,b*y*(1-y))/(b*(a+1))**3
          -mean_h(mp,t,M,a,a,0)/(a*(b+1))**3)


def outer_triangle(mp,t,nu):
    t,nu=map(mp.mpf,(t,nu));A=1+t;root=mp.sqrt(t*t+4*t*nu*(1-nu)+4)
    delta=2*(t*nu*(1-nu)+1)/(root+t)
    return -delta*(2*A+delta)/(A*A*(A+delta)**2)
