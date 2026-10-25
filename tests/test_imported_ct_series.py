import copy,json,unittest
from fractions import Fraction
from unittest.mock import patch
from phi3lab.ct_series_a import build_family
from phi3lab.ct_series_b import ActionLibrary
from phi3lab.ct_series_verify import compare,validate
from phi3lab.ct_series_model import signature
from phi3lab.ct_identity import graph_id
from phi3lab.graphs import Graph


class SeriesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib=ActionLibrary();cls.pairs={}
        for L,E in [(1,2),(1,3),(2,2),(2,3),(3,2),(3,3)]:
            cls.pairs[L,E]=(build_family(L,E),cls.lib.family(L,E))

    def test_complete_low_and_three_loop_families(self):
        for (L,E),(a,b) in self.pairs.items():
            v=compare(a,b)
            self.assertTrue(v['all_equal'],(L,E,v['differences'][:1]))
            self.assertEqual(v['checked_cells'],len(a['rows'])*(len(a['columns'])+1+len(a['overall_columns'])))
            self.assertTrue(all(r['henge_mapping_verified'] and r['forest_oracle_equal'] for r in a['rows']))

    def test_path_b_is_independent_of_all_a_generators(self):
        names=['phi3lab.graphs.generate','phi3lab.graphs.divergent_subgraphs','phi3lab.graphs.contract',
               'phi3lab.higher_graphs.generate','phi3lab.ct_series_a.build_family','phi3lab.ct_series_a.recursive_expansion',
               'phi3lab.renormalization.insertion_columns']
        from contextlib import ExitStack
        with ExitStack() as stack:
            for name in names:stack.enter_context(patch(name,side_effect=AssertionError(name)))
            b=ActionLibrary().family(3,3)
        self.assertTrue(compare(self.pairs[3,3][0],b)['all_equal'])

    def test_two_point_mass_and_kinetic_overall_are_distinct(self):
        a,b=self.pairs[2,2]
        self.assertEqual(set(a['overall_columns']),{'C2mL2','C2zL2'})
        for r in a['rows']:
            mass={t['symbol'] for t in r['terms'] if t['column']=='C2mL2'}
            kinetic={t['symbol'] for t in r['terms'] if t['column']=='C2zL2'}
            self.assertTrue(mass and kinetic);self.assertTrue(mass.isdisjoint(kinetic))
            for t in r['terms']:
                if t['kind']=='overall':
                    self.assertEqual(t['taylor']['degree'],2)
                    self.assertEqual(t['taylor']['operator_degree'],2 if t['column']=='C2zL2' else 0)

    def test_low_counterterms_are_used_recursively(self):
        b=self.pairs[3,3][1]
        events=[e for r in b['rows'] for t in r['terms'] for c in t['contributions'] for e in c.get('gluing_events',[])]
        self.assertTrue(any(Graph.from_dict(e['source']).loops==2 for e in events))
        self.assertTrue(any(len(c.get('gluing_events',[]))>1 for r in b['rows'] for t in r['terms'] for c in t['contributions']))
        self.assertTrue(any(len(e['source_nodes'])>1 for e in events))

    def test_external_host_momenta_and_conservation_are_recorded(self):
        b=self.pairs[3,2][1];found=False
        for r in b['rows']:
            for t in r['terms']:
                for c in t['contributions']:
                    for e in c.get('gluing_events',[]):
                        vectors=list(e['source_boundary_in_host_basis'].values());basis=e['host_routing']['basis']
                        self.assertTrue(all(sum(v[i] for v in vectors)==0 for i in range(len(basis))))
                        found|=any(v[i] for v in vectors for i,k in enumerate(basis) if k.startswith('q'))
        self.assertTrue(found)

    def test_overlapping_products_absent_but_overlapping_subgraphs_retained(self):
        found=False
        for r in self.pairs[3,3][0]['rows']:
            subs=r['recursive_ast']['proper_parts']
            for i,a in enumerate(subs):
                for b in subs[i+1:]:
                    x,y=set(a['edge_ids']),set(b['edge_ids'])
                    if x&y and not(x<=y or y<=x):
                        found=True
                        self.assertFalse(any(x in [set(n['edges']) for n in t['nodes']] and y in [set(n['edges']) for n in t['nodes']] for t in r['terms']))
        self.assertTrue(found)

    def test_coefficient_sign_and_multiplicity_mutations_rejected(self):
        a,b=self.pairs[2,2]
        for factor in (-1,2):
            bad=copy.deepcopy(b);bad['rows'][0]['terms'][0]['coefficient']=str(Fraction(bad['rows'][0]['terms'][0]['coefficient'])*factor)
            with self.assertRaises(ValueError):compare(a,bad)

    def test_label_type_loop_and_mapping_mutations_rejected(self):
        a,b=self.pairs[2,2]
        for mode in ('label','type','loop','momentum','halfedge'):
            bad=json.loads(json.dumps(b))
            if mode=='label':bad['rows'][0]['graph']['external'][0][0]=7
            elif mode=='type':bad['columns'][0]['host']['vertex_kinds'][0]='wrong'
            elif mode=='loop':next(t for r in bad['rows'] for t in r['terms'] if t['nodes'])['nodes'][0]['loops']+=1
            else:
                c=next(c for r in bad['rows'] for t in r['terms'] for c in t['contributions'] if c.get('gluing_events'))
                if mode=='momentum':next(iter(c['gluing_events'][0]['source_boundary_in_host_basis'].values()))[0]+=1
                else:c['canonical_map']['mapping']['half_edges']['0:0'][1]^=1
            with self.assertRaises(ValueError,msg=mode):compare(a,bad)


if __name__=='__main__':unittest.main()
