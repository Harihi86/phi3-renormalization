"""Graph-driven recursive BPH Schwinger-Gaussian compiler.

No topology-specific integration kernel is imported. All original edge labels
remain common to the recursive words. A compatible tree makes each selected
subgraph's loop space explicit; formal boundary scales implement its *total*
external-momentum Taylor polynomial before Gaussian integration. The overall
scale is integrated analytically after the overall zero-momentum subtraction.
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from functools import lru_cache
from itertools import combinations, product
from math import ceil, fsum, isfinite, log10, log1p
from .graphs import Graph, components, degrees, is_1pi
from .kinematics import validate_three_point


def _zero(x):
    return isinstance(x, (int, float, Decimal)) and x == 0


_CANCELLATION_ABSOLUTE_TARGET = 0.
_CANCELLATION_RELATIVE_TARGET = 1e-10
_ROUNDOFF_SAFETY = 64


def _cancellation_risk(value, absolute_components, epsilon=2**-52):
    """A conservative trigger, not a rigorous floating-point error bound."""
    return (not isfinite(float(value)) or not isfinite(float(absolute_components))
            or absolute_components * epsilon * _ROUNDOFF_SAFETY
            > _CANCELLATION_ABSOLUTE_TARGET + abs(value) * _CANCELLATION_RELATIVE_TARGET)


@lru_cache(maxsize=12)
def _jet_layout(orders):
    indices = tuple(product(*(range(n + 1) for n in orders)))
    lookup = {p: i for i, p in enumerate(indices)}
    pairs = []
    for k in indices:
        pairs.append(tuple((lookup[a], lookup[tuple(y-x for x,y in zip(a,k))])
                           for a in indices if all(x <= y for x,y in zip(a,k))))
    return indices, lookup, tuple(pairs)


class _Jet:
    """Truncated multivariate Taylor coefficients, scalar or NumPy batch valued."""
    __slots__ = ('orders', 'c')

    def __init__(self, orders, constant=0, coefficients=None):
        self.orders = tuple(orders)
        self.c = list(coefficients) if coefficients is not None else [constant] + [0]*(_jet_size(orders)-1)

    def coerce(self, other):
        return other if isinstance(other, _Jet) else _Jet(self.orders, other)

    def __add__(self, other):
        other=self.coerce(other)
        return _Jet(self.orders, coefficients=[b if _zero(a) else a if _zero(b) else a+b for a,b in zip(self.c,other.c)])
    __radd__=__add__

    def __neg__(self):
        return _Jet(self.orders, coefficients=[-x for x in self.c])
    def __sub__(self, other): return self + (-self.coerce(other))
    def __rsub__(self, other): return self.coerce(other) + (-self)

    def __mul__(self, other):
        if not isinstance(other, _Jet):
            return _Jet(self.orders, coefficients=[0 if _zero(x) else x*other for x in self.c])
        pairs=_jet_layout(self.orders)[2]
        out=[]
        for terms in pairs:
            value=0
            for a,b in terms:
                if not _zero(self.c[a]) and not _zero(other.c[b]):
                    value=value+self.c[a]*other.c[b]
            out.append(value)
        return _Jet(self.orders, coefficients=out)
    __rmul__=__mul__

    def inverse(self):
        inverse0=1/self.c[0]
        remainder=self*inverse0-1
        term=_Jet(self.orders,1)
        result=term
        for _ in range(sum(self.orders)):
            term=term*(-remainder)
            result=result+term
        return result*inverse0

    def __truediv__(self,other):
        return self*other.inverse() if isinstance(other,_Jet) else _Jet(self.orders,coefficients=[0 if _zero(x) else x/other for x in self.c])

    def __pow__(self, power):
        if power<0:return self.inverse()**(-power)
        result=_Jet(self.orders,1)
        base=self
        while power:
            if power&1:result=result*base
            base=base*base;power//=2
        return result

    def logarithm_shift(self, mass, array_module=None):
        h0=self.c[0]
        logarithm=(array_module.log1p(h0/mass) if array_module is not None else log1p(h0/mass))
        z=(self-h0)/(mass+h0)
        result=_Jet(self.orders,logarithm)
        term=_Jet(self.orders,1)
        for n in range(1,sum(self.orders)+1):
            term=term*z
            result=result+term*((-1)**(n+1))/n
        return result


def _jet_size(orders):return len(_jet_layout(tuple(orders))[0])


def _variable(orders, powers):
    index=_jet_layout(tuple(orders))[1]
    coeff=[0]*len(index)
    if powers in index:coeff[index[powers]]=1
    return _Jet(orders,coefficients=coeff)


def _proper_parts(g):
    """Find embedded divergent edge subgraphs without factorial symmetry work."""
    found=[]
    for size in range(1,len(g.edges)):
        for ids in combinations(range(len(g.edges)),size):
            vertices=sorted({v for i in ids for v in g.edges[i]})
            loops=size-len(vertices)+1
            if loops<1:continue
            omega=6*loops-2*size
            if omega<0:continue
            edges=tuple(g.edges[i] for i in ids)
            if not is_1pi(vertices,edges):continue
            ports=3*len(vertices)-2*size
            found.append({'edge_ids':list(ids),'vertices':vertices,'loops':loops,
                          'external_legs':ports,'omega':omega})
    return found


def _recursive_words(parts):
    """Expand Cγ=-tγ Rbarγ by disjoint spinneys, recursively inner first.

    This directly implements the Bogoliubov recursion; no call to forest
    enumeration or a topology registry occurs here.
    """
    edge_sets=[frozenset(p['edge_ids']) for p in parts]
    vertices=[frozenset(p['vertices']) for p in parts]
    @lru_cache(maxsize=None)
    def expand(available):
        output={frozenset()}
        for size in range(1,len(available)+1):
            for top in combinations(available,size):
                if any(vertices[a]&vertices[b] for a,b in combinations(top,2)):continue
                children=[expand(tuple(b for b in available if edge_sets[b]<edge_sets[a])) for a in top]
                for combination in product(*children):
                    output.add(frozenset(top).union(*combination))
        return tuple(sorted(output,key=lambda word:(len(word),tuple(sorted(word)))))
    return expand(tuple(range(len(parts))))


def _compatible_routing(g, selected):
    edge_order=sorted(range(len(g.edges)),key=lambda i:(min((len(p['edge_ids']) for p in selected if i in p['edge_ids']),default=len(g.edges)+1),i))
    parent=list(range(g.n))
    def root(v):
        while parent[v]!=v:
            parent[v]=parent[parent[v]];v=parent[v]
        return v
    tree=[];chords=[]
    for i in edge_order:
        a,b=g.edges[i];ra,rb=root(a),root(b)
        if ra==rb:chords.append(i)
        else:tree.append(i);parent[ra]=rb
    chords=sorted(chords)
    loops=len(chords);physical=len(g.external)-1;width=loops+physical
    unit=lambda j:tuple(int(k==j) for k in range(width))
    external={label:unit(loops+j) for j,(label,_) in enumerate(g.external[:-1])}
    external[g.external[-1][0]]=tuple(-sum(c[k] for c in external.values()) for k in range(width))
    demand=[[0]*width for _ in range(g.n)]
    for label,v in g.external:
        demand[v]=[a+b for a,b in zip(demand[v],external[label])]
    route={}
    for j,i in enumerate(chords):
        route[i]=unit(j);a,b=g.edges[i]
        demand[a]=[v-c for v,c in zip(demand[a],route[i])]
        demand[b]=[v+c for v,c in zip(demand[b],route[i])]
    remaining=set(tree)
    while remaining:
        degree=[0]*g.n
        for i in remaining:
            a,b=g.edges[i];degree[a]+=1;degree[b]+=1
        leaf=next(v for v,d in enumerate(degree) if d==1)
        i=next(i for i in remaining if leaf in g.edges[i]);a,b=g.edges[i]
        other=b if a==leaf else a
        route[i]=tuple(demand[leaf] if a==leaf else [-x for x in demand[leaf]])
        demand[other]=[x+y for x,y in zip(demand[other],demand[leaf])]
        demand[leaf]=[0]*width;remaining.remove(i)
    if any(any(row) for row in demand):raise ArithmeticError('Routing conservation failed')
    for part in selected:
        inside=set(part['edge_ids']);internal=[j for j,e in enumerate(chords) if e in inside]
        if len(internal)!=part['loops']:raise ArithmeticError('Incompatible subgraph cycle basis')
        if any(route[e][j] for j in internal for e in range(len(g.edges)) if e not in inside):
            raise ArithmeticError('Subgraph loop leaks outside its edges')
    return [list(route[i]) for i in range(len(g.edges))],chords,tree


@dataclass(frozen=True)
class CompiledBPH:
    graph: Graph
    parts: tuple
    words: tuple
    terms: tuple

    def to_dict(self):
        g=self.graph;legs=len(g.external)
        return {'schema':'phi3lab.recursive-gaussian-bph.v1','graph':g.to_dict(),
                'loops':g.loops,'edges':len(g.edges),'external_legs':legs,
                'method':'Graph-derived recursive Bogoliubov BPH; common-edge Schwinger Gaussian Taylor jets',
                'topology_specific_kernel_used':False,'dimension':6,
                'subtraction':'Full zero-momentum Taylor in all independent subgraph boundary momenta',
                'overall_degree':6-2*legs,'proper_subgraphs':list(self.parts),
                'recursive_words':[list(w) for w in self.words],
                'terms':list(self.terms),
                'formula':('0' if legs==1 else 'integral_{alpha>=0,sum alpha=1} sum_w sign(w) sum_{0<=k_gamma<=omega_gamma} [t^k] U_w(t)^(-3) Phi_E(H_w(t),M) d^(I-1)alpha'),
                'gaussian_definition':{'Q':'sum_e alpha_e c_e(t)c_e(t)^T','U':'det Q_loop',
                    'H':'external Gram contraction of Q_external-Q_cross^T Q_loop^-1 Q_cross'},
                'radial_kernel':{'E2':'(M+H) log(1+H/M)-H','E3':'-log(1+H/M)'},
                'normalization':'d^6k/pi^3 per loop; no symmetry weight or coupling',
                'common_regulator':'All original Schwinger alpha_e >= eta until the recursive sum is formed',
                'pointwise_evaluable':True,
                'integration_note':'The common-parameter recursive sum is combined before numerical integration. Numerical QMC uncertainty and UV-face diagnostics are separate from a rigorous total quadrature bound.'}

    def density(self, alpha, invariants=None, m2=1):
        """Evaluate the complete finite projective density at an interior point."""
        if len(alpha)!=len(self.graph.edges) or any(not isfinite(float(a)) or a<=0 for a in alpha):
            raise ValueError('One positive finite parameter is required per original edge')
        total=fsum(alpha)
        if abs(total-1)>1e-10:raise ValueError('Projective parameters must sum to one')
        gram=_external_gram(self.graph,invariants)
        if not isfinite(float(m2)) or m2<=0:raise ValueError('m2 must be positive and finite')
        if self._precision_faces(alpha):return self._checked_density(alpha,invariants,m2)
        condition={}
        try:
            result=self._evaluate(tuple(alpha),gram,m2,None,cancellation_diagnostics=condition)
        except (ArithmeticError,ValueError):
            # Kinematics and mass were checked above; a failed floating
            # Schur/log evaluation is recomputed with the same valid input.
            return self._checked_density(alpha,invariants,m2)
        if (_cancellation_risk(result,condition['absolute_components']) or
                _ROUNDOFF_SAFETY*2**-52*condition.get('gaussian_condition',1)**2>_CANCELLATION_RELATIVE_TARGET):
            return self._checked_density(alpha,invariants,m2)
        return result

    def _checked_density(self,alpha,invariants,mass):
        # The optional executor consumes this exact routing/jet program. The
        # standard-library implementation remains the independent fallback.
        from .automatic_native import try_native_batch_density
        native=try_native_batch_density(self,[alpha],invariants,mass)
        if native is not None:return native[0]
        return float(self._adaptive_decimal(alpha,invariants,mass)[0])

    def density_high_precision(self,alpha,invariants=None,m2=1,*,digits=80,_diagnostics=None):
        """Standard-library Decimal recomputation of the complete recursive sum."""
        if type(digits) is not int or not 30<=digits<=500:raise ValueError('Use30--500 decimal digits')
        if len(alpha)!=len(self.graph.edges) or any(not isfinite(float(a)) or a<=0 for a in alpha):raise ValueError('Positive finite edge parameters are required')
        with localcontext() as context:
            context.prec=digits
            convert=lambda x:Decimal.from_float(float(x))
            aa=tuple(convert(a) for a in alpha)
            gram=tuple(tuple(convert(c) for c in row) for row in _external_gram(self.graph,invariants))
            return self._evaluate(aa,gram,convert(m2),_DecimalMath,cancellation_diagnostics=_diagnostics)

    def _precision_faces(self,alpha):
        all_edges=set(range(len(self.graph.edges)))
        for part in self.parts:
            inside=set(part['edge_ids'])
            if sum(alpha[e] for e in inside)<1e-4 or sum(alpha[e] for e in all_edges-inside)<1e-4:return True
        return False

    def _adaptive_decimal(self,alpha,invariants,mass):
        all_edges=set(range(len(self.graph.edges)));scales=[1.,*alpha]
        for part in self.parts:
            inside=set(part['edge_ids']);scales.extend((sum(alpha[e] for e in inside),sum(alpha[e] for e in all_edges-inside)))
        minimum=min(float(v) for v in scales if v>0)
        digits=min(470,max(80,32+12*ceil(max(0,-log10(minimum)))))
        value=self.density_high_precision(alpha,invariants,mass,digits=digits)
        while True:
            higher=min(500,digits+30)
            details={}
            checked=self.density_high_precision(alpha,invariants,mass,digits=higher,_diagnostics=details)
            with localcontext() as context:
                context.prec=higher
                scale=max(abs(checked),abs(value))
                rounding=Decimal(1024)*Decimal(10)**(-higher+5)*details['absolute_components']
                stable=((checked==value==0 and details['absolute_components']==0) or
                        (scale!=0 and abs(checked-value)<=scale*Decimal('1e-25')
                         and rounding<=abs(checked)*Decimal('1e-25')))
            if stable:return checked,higher
            if higher==500:raise ArithmeticError('Boundary density did not stabilize by500 decimal digits')
            digits=higher;value=checked

    def batch_density(self, alpha, invariants=None, m2=1,*,diagnostics=None):
        """Optional NumPy batch, with shape (samples, original edges)."""
        import numpy as np
        alpha=np.asarray(alpha,dtype=float)
        if alpha.ndim!=2 or alpha.shape[1]!=len(self.graph.edges) or np.any(alpha<=0) or not np.all(np.isfinite(alpha)):
            raise ValueError('Expected positive sample-by-edge parameters')
        if not np.allclose(alpha.sum(axis=1),1,atol=1e-12,rtol=1e-12):raise ValueError('Projective parameters must sum to one')
        if not isfinite(float(m2)) or m2<=0:raise ValueError('m2 must be positive and finite')
        extended=alpha.astype(np.longdouble)
        condition={}
        result=np.asarray(self._evaluate(tuple(extended[:,j] for j in range(alpha.shape[1])),_external_gram(self.graph,invariants),m2,np,cancellation_diagnostics=condition),dtype=np.longdouble)
        absolute=np.asarray(condition['absolute_components'],dtype=np.longdouble)
        cancellation=(~np.isfinite(result)|~np.isfinite(absolute)
                      |(absolute*np.finfo(np.longdouble).eps*_ROUNDOFF_SAFETY
                        >_CANCELLATION_ABSOLUTE_TARGET+np.abs(result)*_CANCELLATION_RELATIVE_TARGET))
        gaussian_condition=np.asarray(condition.get('gaussian_condition',1),dtype=np.longdouble)
        gaussian_risk=(~np.isfinite(gaussian_condition)|
                       (np.finfo(np.longdouble).eps*_ROUNDOFF_SAFETY*gaussian_condition**2
                        >_CANCELLATION_RELATIVE_TARGET))
        flagged=np.zeros(len(alpha),dtype=bool)
        all_edges=set(range(alpha.shape[1]))
        for part in self.parts:
            ids=part['edge_ids'];other=sorted(all_edges-set(ids))
            flagged|=alpha[:,ids].sum(axis=1)<1e-4
            if other:flagged|=alpha[:,other].sum(axis=1)<1e-4
        faces=flagged.copy();flagged|=cancellation|gaussian_risk
        highest=0;decimal_points=0;native_info={}
        positions=np.flatnonzero(flagged)
        if len(positions):
            from .automatic_native import try_native_batch_density
            native=try_native_batch_density(self,alpha[positions].tolist(),invariants,m2,diagnostics=native_info)
            if native is not None:
                result[positions]=native
                highest=native_info.get('maximum_decimal_digits',0)
                decimal_points=native_info.get('decimal_fallback_points',0)
            else:
                for j in positions:
                    checked,digits=self._adaptive_decimal(alpha[j],invariants,m2)
                    result[j]=np.longdouble(str(checked));highest=max(highest,digits)
                decimal_points=len(positions)
        if diagnostics is not None:
            diagnostics['precision_recomputed_points']=diagnostics.get('precision_recomputed_points',0)+len(positions)
            diagnostics['decimal_recomputed_points']=diagnostics.get('decimal_recomputed_points',0)+decimal_points
            diagnostics['maximum_decimal_digits']=max(diagnostics.get('maximum_decimal_digits',0),highest)
            diagnostics['native_recomputed_points']=diagnostics.get('native_recomputed_points',0)+native_info.get('native_checked_points',0)
            if native_info.get('native_checked_points'):
                diagnostics['native_precision_pair']=native_info['native_precision_pair']
                diagnostics['native_library_sha256']=native_info['native_library_sha256']
                diagnostics['native_precision_failed_points']=diagnostics.get('native_precision_failed_points',0)+native_info.get('native_precision_failed_points',0)
                diagnostics['native_double_underflow_points']=diagnostics.get('native_double_underflow_points',0)+native_info.get('native_double_underflow_points',0)
                accepted=[check for check,status in zip(native_info.get('native_relative_checks',[]),native_info.get('native_status',[])) if status==0 and isfinite(check)]
                diagnostics['maximum_native_relative_check']=max(diagnostics.get('maximum_native_relative_check',0),max(accepted,default=0))
            diagnostics['maximum_working_decimal_digits']=max(diagnostics.get('maximum_working_decimal_digits',0),highest,100 if native_info.get('native_checked_points') else 0)
            diagnostics['recomputed_points_precision_rechecked']=True
            diagnostics['base_mantissa_bits']=int(np.finfo(np.longdouble).nmant+1)
            diagnostics['face_threshold']=1e-4
            diagnostics['uv_face_flagged_points']=diagnostics.get('uv_face_flagged_points',0)+int(faces.sum())
            diagnostics['cancellation_flagged_points']=diagnostics.get('cancellation_flagged_points',0)+int(np.count_nonzero(cancellation))
            diagnostics['gaussian_condition_flagged_points']=diagnostics.get('gaussian_condition_flagged_points',0)+int(np.count_nonzero(gaussian_risk))
            diagnostics['maximum_gaussian_condition']=max(diagnostics.get('maximum_gaussian_condition',0),min(1e300,float(np.max(gaussian_condition))))
            diagnostics['gaussian_condition_trigger']='64*eps*max(trace(Q)/pivot, original_external_scale/Schur_H)^2 > relative_target; heuristic'
            diagnostics['cancellation_trigger_absolute_target']=_CANCELLATION_ABSOLUTE_TARGET
            diagnostics['cancellation_trigger_relative_target']=_CANCELLATION_RELATIVE_TARGET
            diagnostics['roundoff_trigger_safety_factor']=_ROUNDOFF_SAFETY
            diagnostics['roundoff_trigger_is_rigorous_error_bound']=False
        return result

    def _evaluate(self,alpha,gram,mass,xp,*,cancellation_diagnostics=None):
        g=self.graph;legs=len(g.external)
        if legs==1:
            if cancellation_diagnostics is not None:
                cancellation_diagnostics['absolute_components']=0
                cancellation_diagnostics['gaussian_condition']=1
            return alpha[0]*0
        results=[]
        absolute_components=0
        gaussian_condition=1
        monitor_gaussian=cancellation_diagnostics is not None and (xp is None or hasattr(xp,'maximum'))
        def update_condition(numerator,denominator):
            nonlocal gaussian_condition
            if xp is None:
                ratio=numerator/denominator if denominator else (1 if not numerator else float('inf'))
                gaussian_condition=max(gaussian_condition,abs(ratio))
            else:
                with xp.errstate(divide='ignore',invalid='ignore'):
                    ratio=xp.where(denominator!=0,numerator/denominator,xp.where(numerator==0,1,float('inf')))
                gaussian_condition=xp.maximum(gaussian_condition,abs(ratio))
        for term in self.terms:
            orders=tuple(term['jet_orders']);n=len(orders);width=g.loops+len(gram)
            matrix=[[_Jet(orders) for _ in range(width)] for _ in range(width)]
            rows=[]
            for row in term['scaled_routing']:
                rows.append([_variable(orders,tuple(v['powers']))*v['coefficient'] if v['coefficient'] else _Jet(orders) for v in row])
            for ae,row in zip(alpha,rows):
                for i in range(width):
                    for j in range(i,width):
                        value=row[i]*row[j]*ae
                        matrix[i][j]=matrix[i][j]+value
            for i in range(width):
                for j in range(i):matrix[i][j]=matrix[j][i]
            if monitor_gaussian:
                gaussian_scale=sum(abs(matrix[i][i].c[0]) for i in range(g.loops))
                external_scale=sum(abs(matrix[g.loops+i][g.loops+j].c[0]*gram[i][j])
                                   for i in range(len(gram)) for j in range(len(gram)))
            determinant=_Jet(orders,1)
            for k in range(g.loops):
                pivot=matrix[k][k]
                if monitor_gaussian:update_condition(gaussian_scale,abs(pivot.c[0]))
                determinant=determinant*pivot
                inverse=pivot.inverse()
                for i in range(k+1,width):
                    factor=matrix[i][k]*inverse
                    for j in range(i,width):
                        matrix[i][j]=matrix[i][j]-factor*matrix[k][j]
                        matrix[j][i]=matrix[i][j]
            h=_Jet(orders)
            for i in range(len(gram)):
                for j in range(len(gram)):
                    h=h+matrix[g.loops+i][g.loops+j]*gram[i][j]
            if monitor_gaussian:update_condition(external_scale,abs(h.c[0]))
            phi=_radial_kernel(h,mass,legs,xp)
            value=determinant**(-3)*phi
            results.append(term['sign']*sum(value.c))
            if cancellation_diagnostics is not None:
                absolute_components=absolute_components+sum(abs(coefficient) for coefficient in value.c)
        if cancellation_diagnostics is not None:
            cancellation_diagnostics['absolute_components']=absolute_components
            cancellation_diagnostics['gaussian_condition']=gaussian_condition
        if xp is None:return fsum(results)
        if not hasattr(xp,'asarray'):return sum(results)
        # Pair recursive terms sample-wise rather than summing integrated pieces.
        return xp.sum(xp.asarray(results),axis=0)


def _external_gram(g,invariants):
    legs=len(g.external)
    if invariants is None:invariants=(1,)*(0 if legs==1 else 1 if legs==2 else 3)
    if legs==1:
        if len(invariants):raise ValueError('A one-point graph has no independent external invariant')
        return ()
    if legs==2:
        if len(invariants)!=1 or not isfinite(float(invariants[0])) or invariants[0]<0:raise ValueError('One nonnegative two-point invariant is required')
        return ((float(invariants[0]),),)
    s1,s2,s3=validate_three_point(invariants)
    cross=fsum((s3/2,-s1/2,-s2/2))
    return ((s1,cross),(cross,s2))


def _radial_kernel(h,mass,legs,xp):
    if legs==3:return -h.logarithm_shift(mass,xp)
    h0=h.c[0];z=h0/mass
    logarithm=xp.log1p(z) if xp is not None else log1p(z)
    if xp is None or not hasattr(xp,'where'):
        if abs(z)<.125:
            if isinstance(z,Decimal):
                term=mass*z*z/2;constant=term;n=2
                while term:
                    n+=1;term=term*(-z)*(n-2)/n
                    next_value=constant+term
                    if next_value==constant:break
                    constant=next_value
            else:
                power=z*z;constant=0
                for n in range(2,30):
                    constant+=mass*((-1)**n)*power/(n*(n-1));power*=z
        else:constant=(mass+h0)*logarithm-h0
    else:
        small=xp.abs(z)<.125
        zz=xp.where(small,z,0);power=zz*zz;series=xp.zeros_like(z)
        for n in range(2,30):
            series+=mass*((-1)**n)*power/(n*(n-1));power*=zz
        constant=xp.where(small,series,(mass+h0)*logarithm-h0)
    delta=h-h0
    result=_Jet(h.orders,constant)+delta*logarithm
    term=delta
    for n in range(2,sum(h.orders)+1):
        term=term*delta
        result=result+term*((-1)**n)/(n*(n-1))*(mass+h0)**(1-n)
    return result


class _DecimalMath:
    @staticmethod
    def log1p(z):
        if abs(z)>=Decimal('.125'):return (1+z).ln()
        term=z;value=z;n=1
        while term:
            n+=1;term=term*(-z)*(n-1)/n
            next_value=value+term
            if next_value==value:return value
            value=next_value
        return value


def compile_bph(g):
    """Compile any generated 1--4-loop E=1,2,3 bare trivalent 1PI graph."""
    if not isinstance(g,Graph):g=Graph.from_dict(g)
    return _compile_bph(g)


@lru_cache(maxsize=1024)
def _compile_bph(g):
    if not 1<=g.loops<=4 or not 1<=len(g.external)<=3 or degrees(g)!=(3,)*g.n or g.vertex_kinds!=('phi3',)*g.n or not is_1pi(range(g.n),g.edges):
        raise ValueError('Compiler supports bare connected trivalent 1PI graphs at 1--4 loops and 1--3 external legs')
    parts=tuple(_proper_parts(g));words=_recursive_words(parts);terms=[]
    for word in words:
        chosen=sorted(word,key=lambda j:(len(parts[j]['edge_ids']),j))
        selected=[parts[j] for j in chosen]
        routing,chords,tree=_compatible_routing(g,selected)
        variables=[j for j in chosen if parts[j]['omega']>0]
        orders=tuple(parts[j]['omega'] for j in variables)
        scaled=[]
        for e,row in enumerate(routing):
            output=[]
            for c,coefficient in enumerate(row):
                powers=[0]*len(variables)
                for j in chosen:
                    if e not in parts[j]['edge_ids']:continue
                    is_internal=c<g.loops and chords[c] in parts[j]['edge_ids']
                    if is_internal:continue
                    if parts[j]['omega']==0:coefficient=0;break
                    powers[variables.index(j)]+=1
                output.append({'coefficient':coefficient,'powers':powers})
            scaled.append(output)
        terms.append({'sign':(-1)**len(word),'subgraphs_inner_first':chosen,
                      'tree_edges':tree,'chord_edges':chords,'unscaled_routing':routing,
                      'jet_subgraphs':variables,'jet_orders':list(orders),'scaled_routing':scaled})
    return CompiledBPH(g,parts,words,tuple(terms))
