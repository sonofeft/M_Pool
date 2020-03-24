import unittest
# import unittest2 as unittest # for versions of python < 2.7

"""
        Method                  Checks that
self.assertEqual(a, b)           a == b   
self.assertNotEqual(a, b)        a != b   
self.assertTrue(x)               bool(x) is True  
self.assertFalse(x)              bool(x) is False     
self.assertIs(a, b)              a is b
self.assertIsNot(a, b)           a is not b
self.assertIsNone(x)             x is None 
self.assertIsNotNone(x)          x is not None 
self.assertIn(a, b)              a in b
self.assertNotIn(a, b)           a not in b
self.assertIsInstance(a, b)      isinstance(a, b)  
self.assertNotIsInstance(a, b)   not isinstance(a, b)  

See:
      https://docs.python.org/2/library/unittest.html
         or
      https://docs.python.org/dev/library/unittest.html
for more assert options
"""

import sys, os

here = os.path.abspath(os.path.dirname(__file__)) # Needed for py.test
up_one = os.path.split( here )[0]  # Needed to find m_pool development version
if here not in sys.path[:2]:
    sys.path.insert(0, here)
if up_one not in sys.path[:2]:
    sys.path.insert(0, up_one)

from m_pool.matrix_pool import MatrixPool, Axis

class MyTest(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.MP = MatrixPool()

        epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 'units':'', 'transform':'log10'})
        pcAxis = Axis({'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':'log10'})
        mrAxis = Axis({'name':'mr', 'valueL':[1,2,3], 'units':'', 'transform':''})
        for A in [epsAxis, pcAxis, mrAxis]:
            self.MP.add_axis( A )
        
        M = self.MP.add_matrix( name='cea_isp', units='sec', axisNameL=['eps','pc','mr'] )
        for eps in epsAxis:
            for pc in pcAxis:
                for mr in mrAxis:
                    M.setByName( pc=pc, eps=eps, mr=mr, val=eps+pc+mr )

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        del( self.MP )

    def test_should_always_pass_cleanly(self):
        """Should always pass cleanly."""
        pass

    def test_MP_existence(self):
        """Check that MP exists"""
        result = self.MP

        # See if the self.MP object exists
        self.assertTrue(result)

    def test_interp(self):
        """Check interpolation"""
        M = self.MP.get_matrix_by_name( 'cea_isp' )
        val = M.interp(order=2, pc=200, eps=20, mr=2.0)
        self.assertAlmostEqual(val, 222.0, places=7)

        val = M.interp(order=2, pc=100, eps=10, mr=1.0)
        self.assertAlmostEqual(val, 111.0, places=7)

        val = M.interp(order=2, pc=150, eps=15, mr=1.5)
        self.assertAlmostEqual(val, 164.69899853053423, places=5)

if __name__ == '__main__':
    # Can test just this file from command prompt
    #  or it can be part of test discovery from nose, unittest, pytest, etc.
    unittest.main()

