import unittest
from fractions import Fraction as Q
from convex_radial_envelope import Node,ConvexEnvelope,positive_weight_error_bound

class EnvelopeTests(unittest.TestCase):
    def test_quadratic_bounds_hold_on_entire_pieces(self):
        envelope=ConvexEnvelope(Node(x,x*x-Q(1,100),x*x+Q(1,100)) for x in (Q(0),Q(1,2),Q(1),Q(3),Q(7)))
        pieces=envelope.pieces();self.assertEqual(pieces[0].left,0);self.assertEqual(pieces[-1].right,7)
        for a,b in zip(pieces,pieces[1:]):self.assertEqual(a.right,b.left)
        for piece in pieces:
            # Upper line minus x² is concave: its minimum is at an endpoint.
            for x in (piece.left,piece.right):self.assertGreaterEqual(piece.upper(x)-x*x,0)
            # x² minus lower line is convex: include its exact stationary point.
            points=[piece.left,piece.right];stationary=piece.lower.slope/2
            if piece.left<stationary<piece.right:points.append(stationary)
            for x in points:self.assertGreaterEqual(x*x-piece.lower(x),0)
            for n in range(11):
                x=piece.left+(piece.right-piece.left)*Q(n,10);self.assertEqual(envelope.bounds(x),(piece.lower(x),piece.upper(x)))
    def test_affine_exact_values_recover_function(self):
        envelope=ConvexEnvelope(Node(x,2*x+1,2*x+1) for x in (Q(0),Q(1),Q(2),Q(4)))
        for n in range(33):
            x=Q(n,8);self.assertEqual(envelope.bounds(x),(2*x+1,2*x+1))
    def test_positive_mass_error_has_exact_integral_upper_bound(self):
        envelope=ConvexEnvelope(Node(x,x*x,x*x) for x in (Q(0),Q(1),Q(2),Q(4)));pieces=envelope.pieces()
        masses=[p.right-p.left for p in pieces];bound=positive_weight_error_bound(pieces,masses)
        actual=sum(((p.upper.slope-p.lower.slope)*(p.right**2-p.left**2)/2+(p.upper.intercept-p.lower.intercept)*(p.right-p.left) for p in pieces),Q())
        self.assertGreaterEqual(bound,actual);self.assertGreater(actual,0)
    def test_invalid_intervals_and_extrapolation_rejected(self):
        for nodes in ([Node(0,0,0)],[Node(1,1,1),Node(1,2,2)],[Node(0,3,3),Node(1,1,1)]):
            with self.assertRaises(ValueError):ConvexEnvelope(nodes)
        envelope=ConvexEnvelope([Node(0,0,0),Node(1,1,1)])
        for x in (-1,2):
            with self.assertRaises(ValueError):envelope.bounds(x)
        with self.assertRaises(ValueError):positive_weight_error_bound(envelope.pieces(),[-1])
        with self.assertRaises(ValueError):positive_weight_error_bound(envelope.pieces(),[])
    def test_nonconvex_node_values_detected(self):
        envelope=ConvexEnvelope([Node(0,0,0),Node(1,2,2),Node(2,3,3)])
        with self.assertRaises(ValueError):envelope.pieces()

if __name__=='__main__':unittest.main(verbosity=2)
