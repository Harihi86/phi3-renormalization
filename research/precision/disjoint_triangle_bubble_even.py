"""Exact q -> -q symmetrization before any numerical angular integration.

The positive mass and complete q-space domain make inversion measure-preserving.
This is an integral identity, not equality of an unsymmetrized point value.
"""
import sympy as s


def verify_identities(coefficients=(1,-1,1)):
    """Check rational parity identities, including signs, without floating point."""
    A,c,d,l,m=s.symbols('A c d l m')
    U=(A+c)*(A+d)+l*m;V=(A+c)*m+(A+d)*l
    determinant=((A+c)**2-l*l)*((A+d)**2-m*m)
    Hp=1/(A*A*(A+c+l)*(A+d+m));Hm=1/(A*A*(A+c-l)*(A+d-m))
    a,b,k=coefficients
    He=a*U/(A*A*determinant);Ho=b*V/(A*A*determinant)
    Hd=k*(V*V-U*(A*(c+d)+c*d+l*m))/(A**4*determinant)
    assert s.cancel(He-(Hp+Hm)/2)==0
    assert s.cancel(Ho-(Hp-Hm)/2)==0
    assert s.cancel(Hd-(Hp+Hm)/2+1/A**4)==0
    M,C,L=s.symbols('M C L')
    delta=(2*M*C+C*C-L*L)/(M*M)
    assert s.cancel((M+C+L)*(M+C-L)/(M*M)-(1+delta))==0
    # The log difference is atanh(L/(M+C)); its derivative and zero value
    # fix the real branch when M+C>|L|. Even log uses the positive product.
    z=s.symbols('z')
    assert s.cancel(s.diff(s.atanh(z),z)-1/(1-z*z))==0 and s.atanh(0)==0
    e,o,f,g=s.symbols('e o f g')
    assert s.expand(((e+o)*(f+g)+(e-o)*(f-g))/2-e*f-o*g)==0
    return {'host_even_odd_rational_identities_exact':True,'host_overall_difference_exact':True,
            'log_product_and_real_atanh_branch_exact':True,'odd_angular_terms_removed_before_evaluation':True,
            'q_inversion_preserves_full_six_dimensional_measure':True,
            'symmetrized_kernel_not_original_pointwise_density':True}


def kernel(ctx,binding,t,qp1,qp2,simplex,gram, *, bubble_value):
    def number(value):
        if isinstance(value,str) and '/' in value:
            num,den=value.split('/');return ctx.mpf(num)/ctx.mpf(den)
        return ctx.mpf(str(value)) if isinstance(value,s.Basic) else ctx.mpf(value)
    t,u,v=map(number,(t,qp1,qp2));x=list(map(number,simplex))
    assert t>=0 and len(x)==3 and all(y>=0 for y in x) and abs(sum(x)-1)<=8*ctx.eps
    g11,g12,g22=map(number,(gram[0][0],gram[0][1],gram[1][1]))
    def ext_norm(a,b):return a*a*g11+2*a*b*g12+b*b*g22
    constants=[];linears=[]
    for a,b in binding['host_external_shifts']:
        a,b=number(a),number(b);constants.append(ext_norm(a,b));linears.append(2*(a*u+b*v))
    A=1+t;c,d=constants;l,m=linears
    assert all(z>0 for z in (A+c+l,A+c-l,A+d+m,A+d-m))
    U=(A+c)*(A+d)+l*m;V=(A+c)*m+(A+d)*l
    determinant=(A+c+l)*(A+c-l)*(A+d+m)*(A+d-m)
    He=U/(A*A*determinant);Ho=-V/(A*A*determinant)
    Hd=(V*V-U*(A*(c+d)+c*d+l*m))/(A**4*determinant)
    V0=ctx.mpf(0);C=ctx.mpf(0);L=ctx.mpf(0)
    for row in binding['triangle_pairs']:
        i,j=row['simplex_pair'];r,a,b=map(number,row['momentum']);weight=x[i]*x[j]
        V0+=weight*r*r*t;C+=weight*ext_norm(a,b);L+=weight*2*r*(a*u+b*v)
    M=1+V0;E=M+C;assert M>0 and E>abs(L)
    even_log_difference=ctx.log1p((2*M*C+C*C-L*L)/(M*M))/2
    odd_log=ctx.atanh(L/E)
    value=-bubble_value*(Hd*ctx.log1p(V0)+He*even_log_difference+Ho*odd_log)
    return {'value':value,'host_even':He,'host_odd':Ho,'host_even_difference':Hd,
            'V0':V0,'constant_external_variance':C,'linear_external_variance':L,
            'even_log_difference':even_log_difference,'odd_log':odd_log}
