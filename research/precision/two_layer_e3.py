"""Graph-derived two-layer symmetric E3 family; no graph-ID lookup or integral values."""
def parameters(tree):
    assert tree['external_legs']==3 and tree['loops']==4
    nested=[c for c in tree['children'] if c['children']]
    leaves=[c for c in tree['children'] if not c['children']]
    assert len(nested)==1
    child=nested[0];assert child['external_legs']==2 and all(not c['children'] and c['loops']==1 and c['external_legs']==2 for c in child['children'])
    assert all(c['loops']==1 and c['external_legs']==2 for c in leaves)
    j=len(leaves);k=len(child['children']);assert (j,k) in ((1,1),(0,2)) and child['loops']==k+1
    assert tree['loops']==j+k+2
    return j,k
