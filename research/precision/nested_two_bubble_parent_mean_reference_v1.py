"""Original routed density, formal jets, axis averages and independent quadratures."""
import mpmath as mp
import nested_two_bubble_reference_v1 as H
import bubble_tail_marginals_v1 as A
M=H.M


def evaluate(point,binding,digits=640,mode='parent-mean',quadrature=False):
 assert binding['id'] in ('L4E3G203','L4E3G64') and mode in ('parent-mean','parent-mean-balanced') and len(point)==5
 ctx=mp.mp.clone();ctx.dps=digits;p=[ctx.mpf(v.numerator)/v.denominator for v in point]
 square_map=lambda u:u*u/(u*u+(1-u)**2)
 map_derivative=lambda u:2*u*(1-u)/(u*u+(1-u)**2)**2
 xi=square_map(p[3]);xi_jac=map_derivative(p[3]);X=1-xi;ell=xi*X
 w=square_map(p[4]) if mode.endswith('balanced') else p[4]
 w_jac=map_derivative(p[4]) if mode.endswith('balanced') else ctx.mpf(1)
 x=X*w;y=xi/(1-x);x0=ctx.sqrt(x)/(ctx.sqrt(x)+ctx.sqrt(1-x))
 s=p[0];scale=1/ctx.sqrt(ell);den=1-s+scale*s;old_s=scale*s/den;radial_jac=scale/den**2
 pp=[old_s,p[1],p[2],x0,y];t=(old_s/(1-old_s))**2;prop=1+t
 factor=radial_jac*xi_jac*X*w_jac/((1-x)*map_derivative(x0))
 full=factor*H.density(ctx,pp,binding,control=False,quadrature=quadrature)
 parent=binding['parent_boundary_shifts'];assert parent[0][0]==-1 and parent[1][0]==0 and parent[0][1:]==parent[1][1:]
 parent3=[*parent,[a-b for a,b in zip(*parent)]];host=binding['host_distinct_routes'];br=binding['disjoint_bubble_boundary']
 norm=lambda rr:ctx.mpf(rr[1]**2+rr[2]**2-rr[1]*rr[2])
 assert norm(parent[0])==1
 weights=[x*(1-x)*y,x*(1-x)*(1-y),(1-x)**2*y*(1-y)]
 V0=ell*t;d=1+V0;K0=M.bubble_K(ctx,x,V0)
 arg=d/x-1;B=H.bubble(ctx,arg);Kp=-x*B/d**2;Kpp=2*x*B/d**3-ctx.diff(lambda u:H.bubble(ctx,u),arg)/d**2
 jfun=lambda u:H.bubble(ctx,u)/(1+u);J0=jfun(t);Jp=ctx.diff(jfun,t);Jpp=ctx.diff(jfun,t,2)
 if quadrature:
  def functions(j,xx):
   if xx==0:return ctx.mpf(0)
   if j==0:return M.bubble_K(ctx,xx,V0)
   arg=d/xx-1;b=H.bubble(ctx,arg);bp=ctx.diff(lambda u:H.bubble(ctx,u),arg)
   kp=-xx*b/d**2;kpp=2*xx*b/d**3-bp/d**2
   return (xx*kp,xx*(1-xx)*kp,xx*xx*kpp)[j-1]
  integrated=[ctx.quad(lambda xx:functions(j,xx),[0,X/2,X]) for j in range(4)]
 else:integrated=A.parent_moments(ctx,xi,t)
 def mul(a,b):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(3)]
 def coefficient(qx,qy,integrated_parent):
  q1,q2=qx,(-qx+ctx.sqrt(3)*qy)/2
  linear=lambda rr:2*rr[0]*(rr[1]*q1+rr[2]*q2)
  jet=[ctx.mpf(1),ctx.mpf(0),ctx.mpf(0)]
  for rr in host:
   l,c=linear(rr)/prop,norm(rr)/prop;jet=mul(jet,[ctx.mpf(1),-l,l*l-c])
  if integrated_parent:
   I0,I1,Ic,I2=integrated;direction=-2*xi*(parent[0][1]*q1+parent[0][2]*q2)
   child=[I0,I1*direction,Ic+I2*direction*direction/2]
  else:
   l=sum(w*linear(rr) for w,rr in zip(weights,parent3));c=sum(w*norm(rr) for w,rr in zip(weights,parent3))
   child=[K0,Kp*l,Kp*c+Kpp*l*l/2]
  jet=mul(jet,child);lb=linear(br)
  jet=mul(jet,[J0,Jp*lb,Jp*norm(br)+Jpp*lb*lb/2])
  return jet[2]
 rad=ctx.sqrt(t*p[1]);theta=p[2]
 point_control=coefficient(rad*ctx.cospi(theta),rad*ctx.sinpi(theta),False)
 root=ctx.sqrt(t);frames=((root,0),(-root,0),(0,root),(0,-root))
 mean=(sum(coefficient(a,b,True) for a,b in frames)+8*coefficient(0,0,True))/12
 radial=2/(s*(1-s))*(t/prop)**3
 return ctx,full-radial*xi_jac*X*w_jac*(1-p[1])*point_control+radial*xi_jac*mean/2
