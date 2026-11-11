"""Independent checks for boundary accounting and the joint optimizer."""
import collections,itertools,random,unittest
from yash_vaccine.constructs.boundaries import junction_peptides,start_peptides,assemble,scan_boundaries
from yash_vaccine.constructs.optimize import solve,objective

class BoundaryTests(unittest.TestCase):
    def test_direct_join_and_full_chain_accounting(self):
        regions=[dict(antigen_id='A',sequence='L'*25),dict(antigen_id='B',sequence='R'*25),dict(antigen_id='C',sequence='A'*25)]
        direct=junction_peptides('L'*25,'','R'*25)
        self.assertEqual(len(direct),34)
        self.assertTrue(all('L' in p and 'R' in p for p in direct))
        links=[dict(name='DIRECT',sequence=''),dict(name='G4S',sequence='GGGGS')]
        sequence,parts=assemble(regions,[0,1,2],links)
        expected=collections.Counter(start_peptides(regions[0]['sequence'])+
              direct+junction_peptides(regions[1]['sequence'],'GGGGS',regions[2]['sequence']))
        actual=collections.Counter(x['peptide'] for x in scan_boundaries(sequence,parts))
        self.assertEqual(actual,expected)
        self.assertEqual([sequence[p['start_aa']:p['end_aa']] for p in parts if p['kind']=='antigen'],[x['sequence'] for x in regions])

    def test_joint_solver_matches_exhaustive_orders_and_linkers(self):
        for seed in [11,23,47]:
            rng=random.Random(seed);n=4
            edges=[dict(left=i,right=j,linker=str(l),sequence='',strong_rank=rng.randrange(3),mild_rank=rng.randrange(6),min_rank=rng.choice([.02,.04,.12,.6,3.1]),added_aa=rng.randrange(6))
                   for i in range(n) for j in range(n) if i!=j for l in range(2)]
            starts=[dict(strong_rank=rng.randrange(2),mild_rank=rng.randrange(3),min_rank=rng.choice([.03,.1,1,4]),added_aa=1) for _ in range(n)]
            lookup={(x['left'],x['right'],int(x['linker'])):x for x in edges}
            brute=min(objective(order,[lookup[(a,b,l)] for a,b,l in zip(order,order[1:],ls)],starts)
                      for order in itertools.permutations(range(n)) for ls in itertools.product(range(2),repeat=n-1))
            order,chosen,meta=solve(n,edges,starts,seconds=5)
            self.assertTrue(meta['all_phases_proven_optimal'])
            self.assertEqual(objective(order,chosen,starts),brute)

if __name__=='__main__':unittest.main()
