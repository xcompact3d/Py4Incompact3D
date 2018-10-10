#!/usr/bin/python3
"""
.. module:: test_tdma
    :synopsis: Unit test for TDMA

.. moduleauthor:: Paul Bartholomew <ptb08@ic.ac.uk>
"""

import unittest
import numpy as np
from Py4Incompact3D.deriv.deriv import tdma

class test_tdma(unittest.TestCase):
    
    def setUp(self):

        self.N = 100
        self.diag = 1.0
        self.coeff = -0.25
        self.f = 1.0

        f = self.f * np.ones(self.N)
        M = self.diag * np.eye(self.N)
        for i in range(self.N):
            if i < (self.N - 1):
                M[i][i + 1] = self.coeff
            if i > 0:
                M[i][i - 1] = self.coeff

        self.x = np.linalg.solve(M, f)

        self.a = self.coeff * np.ones(self.N)
        self.b = self.diag * np.ones(self.N)
        self.c = self.a
        self.rhs = np.array([[self. f * np.ones(self.N)]]) # XXX TDMA expects a 3D array

    def test_solve(self):

        # Compute solution (without overwriting)
        u = tdma(self.a, self.b, self.c, self.rhs, False)
        self.assertTrue(np.allclose(u, self.x))

if __name__ == "__main__":
    unittest.main()

