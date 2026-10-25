"""Independent exact-input cube reflections applied to the full high-precision density."""
import g521_mixture_reference_v1 as M


def evaluate(point,binding,control=True,digits=300,shape=7):
    assert shape==7;total=None
    for mask in range(8):
        p=list(point)
        for i in range(3):
            if mask&(1<<i):p[3+2*i]=1-p[3+2*i]
        ctx,value=M.evaluate(p,binding,control,digits,shape=6)
        total=value/8 if total is None else total+value/8
    return ctx,total
