"""Symmetric E3 one-branch bubble-chain reduction; diagnostic implementation."""
import bubble_chain_candidate as B

def remainder(ctx,t):
    """Positive y-integrated magnitude of the overall zeroth-order remainder.

    The only subtraction has a bounded cancellation ratio <=4/9.
    M=p1^2=p2^2=p3^2=1 is part of this specialization.
    """
    t=ctx.mpf(t);S=ctx.sqrt(t*t+4);q=4/(S+t);D=t+1;H=t+2+S;J=3*S+t+2
    N=8*D*D+4*D*q*(3*D-2)+2*q*q*(8*D-1)+3*q**3
    first=N/(D*D*H*H*J)
    w2=48*t*t/(H*H*J*J);v=w2/(1+w2)
    term=ctx.mpf(1)/3;series=term;n=0
    while True:
        term*=v*(2*n+2)/(2*n+5);series+=term;n+=1
        if not term or abs(term)<ctx.eps*abs(series):break
        if n>1024:raise ArithmeticError('Positive triangle series unresolved')
    second=256*(t*t+t+4)*t/(H**3*J**3)*series/(1+w2)
    assert 0<=second<=ctx.mpf(4)/9*first
    return first-second

def direct_y_remainder(ctx,t,y):
    a=y*(1-y);D=t+1;S=ctx.sqrt(t*t+4+4*t*a)
    delta=2*(1+t*a)/(t+S)
    return delta*(2*D+delta)/(D*D*(D+delta)**2)

def radial(ctx,t,k=3):
    return -t*t/2*B.bubble(ctx,t)**k/(t+1)**(k+1)*remainder(ctx,t)
