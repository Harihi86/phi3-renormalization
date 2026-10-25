"""Full radial-domain diagnostic of the already proved four-dimensional kernel."""
import numpy as np


def bubble(t):
    result=np.empty_like(t);small=t<np.longdouble('.5')
    x=t[small];term=x*x/60;total=term.copy()
    for n in range(2,30):term=-term*x*(n-1)/(2*(2*n+3));total+=term
    result[small]=total
    x=t[~small];result[~small]=(x+4)/3*(np.arcsinh(np.sqrt(x)/2)/np.sqrt(x/(x+4))-1)-x/9
    assert (result>=0).all();return result


def moments(fc,slope,curvature):
    beta=slope/(2*fc);gamma=curvature/(4*fc)
    rho=(np.sqrt(beta*beta+4*gamma)+np.abs(beta))/2
    assert (fc>0).all() and (gamma>=0).all() and (rho<np.longdouble(2)/3).all()
    coefficient=beta*beta+2*gamma;prior=np.full_like(fc,2);current=coefficient.copy();total=np.zeros_like(fc)
    for n in range(1,65):
        total+=current/(2*n*(2*n+1));prior,current=current,coefficient*current-gamma*gamma*prior
    return total


def evaluate(binding,points,chart,radial_power=2):
    assert chart in (0,1) and radial_power in (1,2)
    s,z,theta,u=np.asarray(points,dtype=np.longdouble).T
    t=(s/(1-s))**radial_power
    if chart==0:
        L=np.log1p(t/2);y=np.expm1(u*L)/t;jac=L*np.exp(u*L)/t
    else:y=(1+u)/2;jac=np.full_like(y,np.longdouble(1)/2)
    co,si=np.cos(np.pi*theta),np.sin(np.pi*theta);second=(-co+np.sqrt(np.longdouble(3))*si)/2
    projection={};norms={};inside=None
    for row in binding['triangle_pairs']:
        i,j=row['simplex_pair'];c,a,b=map(int,row['momentum']);norm=a*a-a*b+b*b
        if i==0:projection[j]=c*(a*co+b*second);norms[j]=norm
        else:inside=norm
    assert set(projection)=={1,2} and inside==1
    alpha=y*(1-y);curvature=(1-y)**2;M=1+t*alpha
    Cc=alpha*(norms[1]+norms[2])/2+curvature/4;Cd=alpha*(norms[1]-norms[2])
    q=np.sqrt(t*z);Lc=q*alpha*(projection[1]+projection[2]);Ld=2*q*alpha*(projection[1]-projection[2])
    plus=moments(M+Cc+Lc,Cd+Ld,curvature);minus=moments(M+Cc-Lc,Cd-Ld,curvature)
    c,l=Cc/M,Lc/M;even=np.log1p(2*c+c*c-l*l)/2-(plus+minus)/2
    odd=np.arctanh(Lc/(M+Cc))-(plus-minus)/2
    shifts=binding['host_external_shifts'];h=int(shifts[0][0])*co+int(shifts[0][1])*second;k=int(shifts[1][0])*co+int(shifts[1][1])*second
    A=1+t;X=2+t;eta=4*(t/X)/X*z;beta=(2*A+1)/X**2;da=1-eta*h*h;db=1-eta*k*k
    U=1+eta*h*k;Dnum=eta*(h+k)**2-U*(beta+eta*h*k);AX=A/X
    prefactor=(t/A)**3*(bubble(t)/A)/da/db;odd_coefficient=-AX*AX*(2*q/X)*(h+k)
    # Native log2(t) density times dx/ds; log(2) cancels exactly.
    value=-radial_power/(s*(1-s))*(1-z)*(1-y)*jac*prefactor*(Dnum*np.log(M)+AX*AX*U*even+odd_coefficient*odd)
    assert np.isfinite(value).all()
    return value
