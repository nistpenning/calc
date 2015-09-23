__author__ = 'sbt'

import unittest

from scipy.constants import pi
import numpy as np
import importlib

import mode_analysis_code as mac
importlib.reload(mac)

class TestHexLattice(unittest.TestCase):

    def test_hex_lattice(self):
        # confirm that a lattice with the right numbe of points is generated
        for s in range(0, 10):
            xs, ys = mac.ModeAnalysis.hex_lattice(shells=s)
            npoints = 1 + 6 * np.sum(range(1, s + 1))
            self.assertTrue(len(xs) == npoints)
            self.assertTrue(len(xs) == npoints)


class TestCalculationConsistency(unittest.TestCase):
    def test_073014_axial_frequency_consistency(self):
        """
        Checks to see if the calculation of the axial frequency wz agrees with the
        experimental results from July 30, 2014's data set. This could break if you
        redefined procedures in the Mode Analysis init method.
        """
        print("Beginninng test 07/30/14")
        a = mac.ModeAnalysis(N=91, Vtrap=[-0.0, -203.0, -421.7], Ctrap=1.0, frot=177.0, Vwall=0.10,
                         wall_order=2)
        a.run()
        self.assertTrue(a.wz / (2 * pi) * .90 <= .853E6 <= a.wz / (2 * pi) * 1.10)

        self.assertTrue(
            max(a.axialEvalsE) / (2 * pi) * .90 <= .853E6 <= max(a.axialEvalsE) / (2 * pi) * 1.10)
        print("Passed test 07/30/14 1")
        b = mac.ModeAnalysis(N=91, Vtrap=[-0.0, -423.0, -860.0], Ctrap=1.0, frot=177.0, Vwall=0.10,
                         wall_order=2)
        b.run()
        self.assertTrue(b.wz / (2 * pi) * .90 <= 1.253E6 <= b.wz / (2 * pi) * 1.10)
        self.assertTrue(max(b.axialEvalsE) / (2 * pi) * .90 <= 1.253E6 <= max(b.axialEvalsE) /
                        (2 * pi) * 1.10)
        print("Passed test 07/30/14 2")

        c = mac.ModeAnalysis(N=91, Vtrap=[-0.0, -863.0, -1740.0], Ctrap=1.0, frot=177.0, Vwall=0.10,
                         wall_order=2)
        c.run()

        self.assertTrue(c.wz / (2 * pi) * .90 <= 1.786E6 <= c.wz / (2 * pi) * 1.10)
        self.assertTrue(max(c.axialEvalsE) / (2 * pi) * .90 <= 1.786E6 <= max(c.axialEvalsE) /
          (2 * pi) * 1.10)
        print("Passed test 073014 3")

        d = mac.ModeAnalysis(N=91, Vtrap=[-0.0, -1194.0, -2047.0], Ctrap=1.0, frot=177.0, Vwall=0.10,
                         wall_order=2)
        d.run()
        self.assertTrue(d.wz / (2 * pi) * .90 <= 1.865E6 <= d.wz / (2 * pi) * 1.10)
        self.assertTrue(max(d.axialEvalsE) / (2 * pi) * .90 <= 1.865E6 <= max(d.axialEvalsE) /
          (2 * pi) * 1.10)
        print("Passed test 073014 4")

    def test_043015_axial_mode_consistency(self):
        """
        Checks to see if the axial frequency (which is identically equal to the COM mode eigenfrequency)
        agrees with the values obtained for varying frot in the lab on April 30th, 2015's data set.
        """
        print("beginning test 043015")
        a = mac.ModeAnalysis(N=331, Vtrap=[-0.0, -1750.0, -2000.0], Ctrap=1.0, frot=177.5, Vwall=2.13,
                         wall_order=2)
        a.run()
        # Using units of hertz
        self.assertTrue(a.wz / (2 * pi) * .90 <= 1.576E6 <= 1.10 * a.wz / (2 * pi))
        self.assertTrue(
            max(a.axialEvalsE) / (2 * pi) * .90 <= 1.576E6 <= 1.10 * max(a.axialEvalsE) / (2 * pi))
        print("Passed test 043015")
        # b = ModeAnalysis(N=331 ,Vtrap=[-0.0,-1750.0,-2000.0], Ctrap = 1.0, frot=180.0,
        # Vwall= 2.13, wall_order=2)
        # b.run_quiet()
        # self.assertTrue(b.wz/(2*pi)*.90<= 1.578E6 and 1.578E6 <= 1.10*b.wz/(2*pi))

        # def test_051415_plane_stability(self):
        # Checks to make sure that we start to see a 1-2 plane transistion
        #  at appropriate frequencies.

        #   for w in np.linspace(187.5, 210, 20):
        #      a = ModeAnalysis(N=331, Vtrap=[0.0,-1750.0,-2000.0], Ctrap = 1.0,
        # frot=w, Vwall= 2.13, wall_order=2)
        #     a.run()


if __name__ == '__main__':
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculationConsistency)
    #unittest.TextTestRunner(verbosity=1).run(suite)
    unittest.main(exit=False)
