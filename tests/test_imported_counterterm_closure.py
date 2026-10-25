"""Meaningful exact closure tests; no large numerical integrals."""
import copy
import json
from fractions import Fraction as F
from itertools import permutations
from math import prod
from unittest.mock import patch
import unittest
from phi3lab.graphs import Graph, generate, canonical, symmetry, is_1pi, divergent_subgraphs
from phi3lab.higher_graphs import canonical as fast_canonical, symmetry as fast_symmetry
from phi3lab.ct_graph import identity, isomorphism, momentum_transport
from phi3lab.ct_wick import enumerate_wick, action_sectors
from phi3lab.ct_bph import build_path_a, recursion, henge_trace, taylor_trace
from phi3lab.ct_action import build_path_b, glue
from phi3lab.ct_verify import compare_paths, check_henge, bridge_oracle
from phi3lab.examples import taylor_product, bubble_taylor, delta, dot


class CountertermClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a=build_path_a(); cls.b=build_path_b()

    def test_all_cells_equal_exactly_including_zeros(self):
        r=compare_paths(self.a,self.b)
        self.assertTrue(r['all_equal']);self.assertFalse(r['numeric_tolerance_used'])
        self.assertEqual(r['checked_cells'],len(self.a['rows'])*(len(self.a['columns'])+2))
        self.assertGreater(r['checked_cells'],r['nonzero_cells'])

    def test_b_cannot_call_parent_or_subgraph_generators(self):
        with (patch('phi3lab.graphs.generate',side_effect=AssertionError('bare catalogue')),
             patch('phi3lab.graphs.divergent_subgraphs',side_effect=AssertionError('A subgraphs')),
             patch('phi3lab.graphs.contract',side_effect=AssertionError('A contraction')),
             patch('phi3lab.renormalization.insertion_columns',side_effect=AssertionError('old columns'))):
            self.assertTrue(compare_paths(self.a,build_path_b())['all_equal'])

    def test_no_target_coefficients_in_path_b_inputs(self):
        import inspect
        self.assertEqual(list(inspect.signature(build_path_b).parameters),[])
        self.assertNotIn('ct_bph',build_path_b.__globals__)

    def test_sign_mutation_rejected(self):
        b=copy.deepcopy(self.b); b['rows'][0]['terms'][0]['coefficient']=str(-F(b['rows'][0]['terms'][0]['coefficient']))
        self.assertFalse(compare_paths(self.a,b)['all_equal'])

    def test_insertion_multiplicity_mutation_rejected(self):
        b=copy.deepcopy(self.b)
        t=next(t for r in b['rows'] for t in r['terms'] if t['kind']=='insertion')
        t['coefficient']=str(2*F(t['coefficient']))
        self.assertFalse(compare_paths(self.a,b)['all_equal'])

    def test_label_mutation_rejected(self):
        b=copy.deepcopy(self.b);b['rows'][0]['graph']['external'][0][0]=4
        with self.assertRaises(ValueError):compare_paths(self.a,b)

    def test_type_mutation_rejected(self):
        b=copy.deepcopy(self.b);b['columns'][0]['host']['vertex_kinds'][0]='wrong_kind'
        with self.assertRaises(ValueError):compare_paths(self.a,b)

    def test_half_edge_mapping_mutation_rejected(self):
        b=copy.deepcopy(self.b);t=next(t for r in b['rows'] for t in r['terms'] if t['kind']=='insertion')
        t['gluing']['expanded_to_canonical']['half_edges']['0:0'][1]^=1
        with self.assertRaises(ValueError):compare_paths(self.a,b)

    def test_momentum_mapping_mutation_rejected(self):
        b=copy.deepcopy(self.b);t=next(t for r in b['rows'] for t in r['terms'] if t['kind']=='insertion')
        t['gluing']['loop_momentum_change']['matrix'][0][0]+=1
        with self.assertRaises(ValueError):compare_paths(self.a,b)

    def test_one_loop_hand_coefficients(self):
        for e,expected in [(1,F(1,2)),(2,F(1,2)),(3,F(1))]:
            d=enumerate_wick(('phi3',)*e,tuple(range(1,e+1)))
            accepted=[r for r in d['records'] if r['accepted']]
            self.assertEqual(len(accepted),1);self.assertEqual(F(accepted[0]['positive_weight']),expected)

    def test_all_low_loop_wick_sets_and_automorphisms(self):
        for loops in (1,2):
            for e in (1,2,3):
                d=enumerate_wick(('phi3',)*(2*loops+e-2),tuple(range(1,e+1)))
                actual={r['id']:F(r['positive_weight']) for r in d['records'] if r['accepted']}
                expected={identity(g):F(1,symmetry(g)) for g in generate(loops,e)}
                self.assertEqual(actual,expected)

    def test_complete_generating_function_pairing_total(self):
        for n,e in [(1,1),(2,2),(3,3),(5,3)]:
            d=enumerate_wick(('phi3',)*n,tuple(range(1,e+1)),True)
            self.assertEqual(d['full_wick_multiplicity'],prod(range(1,3*n+e,2)))
            for r in d['records']:
                self.assertEqual(F(r['positive_weight']),F(1,symmetry(Graph.from_dict(r['graph']))))

    def test_raw_wick_half_edges_are_explicitly_bijective(self):
        d=enumerate_wick(('phi3',)*5)
        for r in d['raw_terms']:
            m=r['witness_stub_to_canonical_halfedge']
            self.assertEqual(len(m),12)
            self.assertEqual(len({tuple(v) for v in m.values()}),12)
            self.assertEqual(set(m),{f'{v}:{s}' for v,s in r['remaining_half_edges']})
            self.assertEqual({tuple(v) for v in m.values()},{(i,end) for i in range(6) for end in (0,1)})

    def test_amputation_does_not_imply_one_pi(self):
        d=enumerate_wick(('phi3',)*5,include_source_pairs=True)
        self.assertTrue(any(r['connected'] and not r['one_pi'] for r in d['records']))
        self.assertTrue(any(r['vacuum_components'] for r in d['records']))
        self.assertTrue(any(r['source_pairs'] for r in d['records']))
        # A tadpole exists even after setting J=0 after one derivative.
        d1=enumerate_wick(('phi3',),(1,),True)
        self.assertTrue(any(r['accepted'] and r['self_loops'] for r in d1['records']))

    def test_all_vertex_edge_endpoint_relabellings(self):
        for g in generate(2,3):
            target=identity(g);fast=fast_canonical(g)[0]
            for p in permutations(range(g.n)):
                h=Graph(g.n,tuple((p[b],p[a]) for a,b in reversed(g.edges)),tuple((l,p[v]) for l,v in g.external))
                self.assertEqual(identity(h),target)
                self.assertEqual(fast_canonical(h)[0],fast)
                self.assertEqual(fast_symmetry(h),symmetry(g))

    def test_colored_canonicalization_and_fixed_labels(self):
        for c in self.b['columns']:
            g=Graph.from_dict(c['host']);self.assertEqual(symmetry(g),fast_symmetry(g))
            for p in permutations(range(g.n)):
                kinds=['']*g.n
                for v,k in enumerate(g.vertex_kinds):kinds[p[v]]=k
                h=Graph(g.n,tuple((p[a],p[b]) for a,b in g.edges),tuple((l,p[v]) for l,v in g.external),tuple(kinds))
                self.assertEqual(identity(h),identity(g));self.assertEqual(fast_canonical(h)[0],fast_canonical(g)[0])
        marked=Graph(2,((0,1),(0,1)),((1,0),(2,1)),('phi3','C3'))
        swapped=Graph(2,marked.edges,((2,0),(1,1)),marked.kinds)
        self.assertNotEqual(identity(marked),identity(swapped))

    def test_external_permutation_orbits_and_weights(self):
        gs=generate(2,3); groups={identity(g,False):g for g in gs}
        self.assertEqual(sum((F(1,symmetry(g)) for g in gs),F()),6*sum((F(1,symmetry(g,False)) for g in groups.values()),F()))
        for r in self.a['rows']:self.assertTrue(set(r['external_permutation_orbit']) <= {identity(g) for g in gs})

    def test_parallel_loop_and_bridge_identities(self):
        self.assertEqual(bridge_oracle(range(2),[(0,1),(0,1)]),[])
        self.assertEqual(bridge_oracle(range(3),[(0,0),(0,1),(1,2),(1,2)]),[1])
        self.assertTrue(is_1pi(range(1),[(0,0)]))
        self.assertFalse(is_1pi(range(2),[(0,0),(0,1)]))
        self.assertEqual(symmetry(Graph(1,((0,0),),((1,0),))),2)

    def test_henge_every_edge_and_convergent_blocks(self):
        has_convergent=False
        for loops in (1,2):
            for legs in (1,2,3):
                for g in generate(loops,legs):
                    records=henge_trace(g)
                    self.assertEqual(len(records),len(g.edges));self.assertTrue(check_henge(g,records))
                    self.assertTrue(all(h['skeleton_loops']==1 for h in records))
                    has_convergent |= any(b['loops'] and not b['needs_uv_subtraction'] for h in records for b in h['blocks'])
        self.assertTrue(has_convergent)

    def test_actual_recursive_ast_and_executable_compiler(self):
        for r in self.a['rows']:
            self.assertTrue(r['checks']['new_recursive_ast_equals_forest_oracle'])
            self.assertTrue(r['checks']['all_spinney_coverage_ok'])
            node=r['recursive_ast'];program=r['executable_taylor_program']
            self.assertEqual(len(node['rbar']),len(program['recursive_words']))
            self.assertEqual([p['edge_ids'] for p in node['proper_parts']],[p['edge_ids'] for p in program['proper_subgraphs']])

    def test_four_point_subdivergence_has_no_overall_counterterm(self):
        # Box with one edge replaced by a one-loop two-point insertion.
        g=Graph(6,((0,4),(4,5),(4,5),(5,1),(1,2),(2,3),(0,3)),((1,0),(2,1),(3,2),(4,3)))
        n=recursion(g);self.assertEqual(n['omega'],-2);self.assertIsNone(n['counterterm'])
        self.assertTrue(n['proper_parts']);self.assertEqual(n['renormalized']['operation'],'identity')

    def test_outer_loop_is_retained_in_subgraph_boundary(self):
        traces=[t['taylor'] for r in self.a['rows'] for t in r['terms'] if t['kind']=='insertion']
        self.assertTrue(all(t['outside_loop_as_boundary'] for t in traces))
        for t in traces:
            self.assertEqual(len(t['inner_loop_basis']),1);self.assertEqual(len(t['outer_loop_basis']),1)
            self.assertEqual(t['degree'],2 if len(t['boundary_ports'])==2 else 0)
            self.assertTrue(t['outside_propagators_retained'])

    def test_exact_total_degree_two_taylor_not_quartic(self):
        q=(F(1),F(2),0,0,0,0);p=(F(1,3),0,F(2,5),0,0,0);z=(0,)*6
        co=taylor_product([(q,z),(q,p)],2)
        self.assertEqual(len(co),3);self.assertEqual(sum(co),bubble_taylor(q,p))
        self.assertEqual(co[2],-dot(p,p)*delta(q)**3+4*dot(q,p)**2*delta(q)**4)
        self.assertEqual(taylor_product([(q,z),(q,tuple(2*x for x in p))],2)[2],4*co[2])

    def test_kinetic_tensor_projection_in_six_dimensions(self):
        q=(F(1),F(2),F(-1),0,0,0);z=(0,)*6;dimension=6
        summed=F()
        for j in range(dimension):
            p=tuple(F(i==j) for i in range(dimension))
            summed+=taylor_product([(q,z),(q,p)],2)[2]
        self.assertEqual(summed/dimension,-delta(q)**3+F(4,dimension)*dot(q,q)*delta(q)**4)

    def test_odd_two_point_term_cancels_by_exact_inversion(self):
        q=(F(1),F(2),0,0,0,0);p=(F(2),F(-1,3),0,0,0,0);z=(0,)*6
        a=taylor_product([(q,z),(q,p)],2)
        b=taylor_product([(tuple(-x for x in q),z),(tuple(-x for x in q),p)],2)
        self.assertEqual(a[1]+b[1],0);self.assertEqual(a[2],b[2])

    def test_tadpole_action_sectors_excluded_by_structure(self):
        sectors=[s for s in self.b['sectors'] if 'C1' in s['action_monomial']]
        self.assertTrue(sectors);self.assertTrue(all(not s['accepted_graphs'] for s in sectors))

    def test_new_local_counterterm_has_unit_wick_host_weight(self):
        self.assertEqual(F(self.b['local_host']['positive_weight']),1)
        for r in self.b['rows']:
            low={}
            for t in r['terms']:
                if t['kind']!='overall':low[t['operator']]=low.get(t['operator'],F())+F(t['coefficient'])
            overall={}
            for t in r['terms']:
                if t['kind']=='overall':overall[t['projection_of']]=overall.get(t['projection_of'],F())+F(t['coefficient'])
            self.assertEqual(overall,{k:-v for k,v in low.items()})


if __name__=='__main__':unittest.main()
