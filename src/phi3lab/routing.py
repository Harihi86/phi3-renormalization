"""Integer momentum routings from a spanning tree and Kirchhoff conservation."""
from fractions import Fraction
from .graphs import connected


def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)


def momentum_routing(g):
    if not connected(range(g.n),g.edges):raise ValueError("Routing requires a connected graph")
    parent=list(range(g.n))
    def root(x):
        while parent[x]!=x:x=parent[x]
        return x
    tree=[];chords=[]
    for i,(a,b) in enumerate(g.edges):
        ra,rb=root(a),root(b)
        if ra==rb:chords.append(i)
        else:tree.append(i);parent[ra]=rb
    labels=[label for label,_ in g.external]
    names=[f"q{i+1}" for i in range(len(chords))]+[f"p{i}" for i in labels[:-1]]
    zero=(0,)*len(names)
    unit=lambda j:tuple(int(i==j) for i in range(len(names)))
    external={label:unit(len(chords)+j) for j,label in enumerate(labels[:-1])}
    total=zero
    for v in external.values():total=add(total,v)
    if labels:external[labels[-1]]=neg(total)
    demand={v:zero for v in range(g.n)}
    for label,v in g.external:demand[v]=add(demand[v],external[label])
    route={}
    for j,i in enumerate(chords):
        route[i]=unit(j);a,b=g.edges[i]
        demand[a]=add(demand[a],neg(route[i]));demand[b]=add(demand[b],route[i])
    remaining=set(tree)
    while remaining:
        degree={v:0 for v in range(g.n)}
        for i in remaining:
            for v in g.edges[i]:degree[v]+=1
        leaf=min(v for v,d in degree.items() if d==1)
        edge=next(i for i in remaining if leaf in g.edges[i]);a,b=g.edges[edge]
        other=b if a==leaf else a
        route[edge]=demand[leaf] if a==leaf else neg(demand[leaf])
        demand[other]=add(demand[other],demand[leaf]);demand[leaf]=zero
        remaining.remove(edge)
    if any(any(v) for v in demand.values()):raise AssertionError("Nonzero residual momentum")
    residual={v:zero for v in range(g.n)}
    for i,(a,b) in enumerate(g.edges):
        residual[a]=add(residual[a],route[i]);residual[b]=add(residual[b],neg(route[i]))
    for label,v in g.external:residual[v]=add(residual[v],neg(external[label]))
    return {"basis":names,"chord_edges":chords,"tree_edges":tree,
            "edge_momenta":{str(i):list(route[i]) for i in range(len(g.edges))},
            "external_momenta":{str(k):list(v) for k,v in external.items()},
            "conservation_ok":all(not any(v) for v in residual.values())}


def momentum_text(coefficients,names):
    parts=[]
    for c,name in zip(coefficients,names):
        if not c:continue
        term=name if abs(c)==1 else f"{abs(c)}*{name}"
        parts.append(("-" if c<0 else "+")+term)
    return "".join(parts).lstrip("+") or "0"
