import unittest
from expression_weights import bounded_fraction,expression_proxy,fusion_proxies
class WeightTests(unittest.TestCase):
 def test_informative_fraction_not_fusion_to_unfused_ratio(self):
  self.assertAlmostEqual(bounded_fraction(80,20),.8)
  self.assertAlmostEqual(expression_proxy(100,.8),80)
  with self.assertRaises(ValueError):expression_proxy(100,4) # FAR is not a fraction.
 def test_missing_and_zero(self):
  self.assertIsNone(bounded_fraction(0,0));self.assertIsNone(expression_proxy(None,.5))
  self.assertEqual(expression_proxy(0,.5),0);self.assertEqual(bounded_fraction(0,20),0)
 def test_two_partners_remain_separate(self):
  d=fusion_proxies(100,200,20,80,20)
  self.assertEqual(d['left_expression_proxy'],20);self.assertEqual(d['right_expression_proxy'],100)
if __name__=='__main__':unittest.main()
