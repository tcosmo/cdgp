import unittest
from cdgp import *

class ExampleTest(unittest.TestCase):
    def test_square(self):
        self.assertEqual(square(5),23) 
