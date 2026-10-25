"""Unified structural identity facade; CT1 identities are never reassigned."""
import hashlib,json
from collections import Counter,defaultdict
from .graphs import Graph
from .ct_graph import identity as ct1_identity, isomorphism as ct1_isomorphism
from .higher_graphs import canonical,symmetry


def graph_id(g,fixed_labels=True):
    # CT1 was bounded to <=5 vertices. Preserve all of that namespace exactly.
    if g.n<=5:return ct1_identity(g,fixed_labels)
    key,_=canonical(g,fixed_labels)
    return 'g2_'+hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()


def mapping(g,target):
    if g.n<=5:return ct1_isomorphism(g,target)
    if canonical(g)[0]!=canonical(target)[0]:raise ValueError('Fixed-label typed graphs are not isomorphic')
    a=Counter(g.edges);b=Counter(target.edges)
    def color(h,v):
        return h.vertex_kinds[v],tuple(l for l,w in h.external if w==v),sum(v in e for e in h.edges),sum(x==y==v for x,y in h.edges)
    pools={v:[w for w in range(target.n) if color(g,v)==color(target,w)] for v in range(g.n)}
    chosen={};used=set()
    def search():
        if len(chosen)==g.n:return True
        v=min((v for v in range(g.n) if v not in chosen),key=lambda v:(sum(w not in used for w in pools[v]),-sum(a[tuple(sorted((v,u)))] for u in chosen),v))
        for w in pools[v]:
            if w in used or any(a[tuple(sorted((v,u)))]!=b[tuple(sorted((w,z)))] for u,z in chosen.items()):continue
            chosen[v]=w;used.add(w)
            if search():return True
            used.remove(w);del chosen[v]
        return False
    if not search():raise ValueError('No adjacency-preserving mapping')
    bins=defaultdict(list)
    for i,e in enumerate(target.edges):bins[e].append(i)
    em={};hm={}
    for i,(v,w) in enumerate(g.edges):
        u,z=chosen[v],chosen[w];j=bins[tuple(sorted((u,z)))].pop(0);flip=int(u>z)
        em[str(i)]=j;hm[f'{i}:0']=[j,flip];hm[f'{i}:1']=[j,1-flip]
    return {'vertices':[chosen[v] for v in range(g.n)],'edges':em,'half_edges':hm,
            'external':{str(l):[l,chosen[v]] for l,v in g.external}}


def record(g,legacy_id=None):
    key,cg=canonical(g)
    return {'id':graph_id(g),'identity_schema':'ct1-sha256-prefix20' if g.n<=5 else 'typed-fixed-label-ir-sha256-v2',
            'canonical_key':key,'graph':g.to_dict(),'canonical_graph':cg.to_dict(),
            'to_canonical':mapping(g,cg),'symmetry':symmetry(g),'legacy_id':legacy_id,
            'fixed_external_labels':True,'coordinates_affect_identity':False,
            'migration':'CT1 ID retained' if g.n<=5 else 'new namespace; no prior CT1 identity for >5 vertices'}
