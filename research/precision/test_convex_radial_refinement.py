import unittest
from fractions import Fraction as Q
from convex_radial_envelope import Node,ConvexEnvelope
from convex_radial_refinement import inspect_envelope,dyadic_input

class RefinementTests(unittest.TestCase):
    def test_affine_has_zero_relative_gap(self):
        env=ConvexEnvelope(Node(x,3*x+1,3*x+1) for x in (0,1,2,4))
        result=inspect_envelope(env,1,4,Q(1,2000));self.assertTrue(result['all_passed']);self.assertEqual(result['maximum_relative_gap'],0)
    def test_exact_quadratic_adaptive_envelope(self):
        xs=[Q(0),Q(1),Q(2),Q(4)]
        for _ in range(12):
            env=ConvexEnvelope(Node(x,x*x,x*x) for x in xs);r=inspect_envelope(env,1,4,Q(1,2000))
            for p in env.pieces():
                if p.right<=1:continue
                a=max(Q(1),p.left);b=p.right
                # For this quadratic, check the exact minimum of x²-lower.
                station=p.lower.slope/2;points=[a,b]+([station] if a<station<b else [])
                for x in points:self.assertGreaterEqual(x*x,p.lower(x))
                for x in (a,b):self.assertGreaterEqual(p.upper(x),x*x)
            if r['all_passed']:break
            xs=sorted(xs+[(xs[i]+xs[i+1])/2 for i in r['failing_intervals']])
        self.assertTrue(r['all_passed']);self.assertLessEqual(r['maximum_relative_gap'],Q(1,2000))
    def test_domain_and_zero_lower_rejected(self):
        env=ConvexEnvelope(Node(x,x*x,x*x) for x in (0,1,2,4))
        for a,b,t in ((0,4,Q(1,10)),(-1,4,Q(1,10)),(1,5,Q(1,10)),(1,4,0)):
            with self.assertRaises(ValueError):inspect_envelope(env,a,b,t)
    def test_exact_native_dyadic_serialization(self):
        for x in (Q(1,1024),Q(3,4096),Q(5,2),Q(2)**49,3*Q(2)**48):
            n,e=dyadic_input(x);self.assertEqual(n*Q(2)**e,x)
        for x in (0,-1,Q(1,3),Q(2)**61):
            with self.assertRaises(ValueError):dyadic_input(x)

if __name__=='__main__':unittest.main(verbosity=2)
