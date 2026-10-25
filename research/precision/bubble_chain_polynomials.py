"""Independent graph-factorized U/F coefficients; no routing determinants."""
from fractions import Fraction as Q
from collections import defaultdict
from itertools import product


def factors(plan, state, invariant):
    pairs=plan['bubble_edges']; chain=plan['chain_edges']; outer=plan['outer_edge']
    variables=state['variables']; U=defaultdict(Q); F=defaultdict(Q)
    for choice in product(*pairs):
        mask=sum(1<<e for e in choice)
        for e in chain+[outer]: U[(mask|(1<<e),0)] += 1
        for e in chain: F[(mask|(1<<outer)|(1<<e),0)] += Q(invariant)
    for j,pair in enumerate(pairs):
        slot=1<<variables.index(j) if j in variables else 0
        for choice in product(*(p for i,p in enumerate(pairs) if i!=j)):
            mask=sum(1<<e for e in (*pair,*choice))
            U[(mask,slot)] += 1
            F[(mask|(1<<outer),slot)] += Q(invariant)
    return dict(U),dict(F)


def original(state, family, invariants):
    out=defaultdict(Q)
    for row in state[family]:
        mask,slot,channel=row[:3]
        c=Q(row[3] if len(row)==4 else 1)
        if family=='F':c*=Q(invariants[channel])
        out[(mask,slot)] += c
    return {key:value for key,value in out.items() if value}


def verify(plan, positive_descriptor, invariants):
    assert len(invariants)==1
    witnesses=[]
    assert len(positive_descriptor['states'])==1<<plan['bubble_count']
    for state in positive_descriptor['states']:
        expected=factors(plan,state,invariants[0])
        actual=tuple(original(state,f,invariants) for f in ('U','F'))
        assert expected==actual,('Independent full original U/F coefficient mismatch',state['word'])
        assert state['sign']==(-1)**len(state['word'])
        witnesses.append({'word':state['word'],'sign':state['sign'],'variables':state['variables'],
            'U':[[mask,slot,str(c)] for (mask,slot),c in sorted(expected[0].items())],
            'F':[[mask,slot,str(c)] for (mask,slot),c in sorted(expected[1].items())]})
    assert {tuple(w['word']) for w in witnesses}=={tuple(w['word']) for w in plan['states']}
    return witnesses
