"""Compress a rigorous affine envelope into certified affine surrogate pieces."""
from fractions import Fraction as Q

def clip(envelope,left,right):
    pieces=[]
    for p in envelope:
        a=max(Q(left),p['left']);b=min(Q(right),p['right'])
        if b<=a:continue
        pieces.append({'left':a,'right':b,'lower':p['lower'],'upper':p['upper']})
    assert pieces and pieces[0]['left']==left and pieces[-1]['right']==right
    assert all(a['right']==b['left'] for a,b in zip(pieces,pieces[1:]))
    return pieces

def midpoint(piece,x):
    return sum(slope*x+intercept for slope,intercept in (piece['lower'],piece['upper']))/2

def certify(source,i,j,epsilon):
    first,last=source[i],source[j];left,right=first['left'],last['right']
    ml,mr=midpoint(first,left),midpoint(last,right);assert ml>0 and mr>0
    slope=(mr-ml)/(right-left);worst=Q()
    for p in source[i:j+1]:
        for x in (p['left'],p['right']):
            lower=p['lower'][0]*x+p['lower'][1];upper=p['upper'][0]*x+p['upper'][1]
            value=ml+slope*(x-left);assert value>0 and 0<lower<=upper
            worst=max(worst,abs(lower-value)/value,abs(upper-value)/value)
            if worst>epsilon:return None
    return {'left':left,'right':right,'mid_left':ml,'mid_right':mr,'source_begin':i,'source_end':j,'maximum_relative_enclosure_error':worst}

def compress(source,epsilon):
    epsilon=Q(epsilon);assert 0<epsilon<1;result=[];i=0
    while i<len(source):
        best=certify(source,i,i,epsilon);assert best is not None
        j=i+1
        while j<len(source):
            candidate=certify(source,i,j,epsilon)
            if candidate is None:break
            best=candidate;j+=1
        result.append(best);i=best['source_end']+1
    assert all(a['right']==b['left'] for a,b in zip(result,result[1:]))
    return result
