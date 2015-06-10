__author__ = 'sbt'

import unittest

from scipy.constants import pi

from crystal_mode_code.mode_analysis_code import ModeAnalysis


class TestCalculationConsistency(unittest.TestCase):

    def test_073014_axial_frequency_consistency(self):
        """
        Checks to see if the calculation of the axial frequency wz agrees with the
        experimental results from July 30, 2014's data set. This could break if you
        redefined procedures in the Mode Analysis init method.
        """
        a = ModeAnalysis(shells=5, Vtrap=[-0.0, -203.0, -421.7], Ctrap=1.0, frot=177.0, Vwall=0.10, wall_order=2)
        self.assertTrue(a.wz / (2 * pi) * .90 <= .853E6 <= a.wz / (2 * pi) * 1.10)

        b = ModeAnalysis(shells=5, Vtrap=[-0.0, -423.0, -860.0], Ctrap=1.0, frot=177.0, Vwall=0.10, wall_order=2)
        self.assertTrue(b.wz / (2 * pi) * .90 <= 1.253E6 <= b.wz / (2 * pi) * 1.10)

        c = ModeAnalysis(shells=5, Vtrap=[-0.0, -863.0, -1740.0], Ctrap=1.0, frot=177.0, Vwall=0.10, wall_order=2)
        self.assertTrue(c.wz / (2 * pi) * .90 <= 1.786E6 <= c.wz / (2 * pi) * 1.10)

        d = ModeAnalysis(shells=5, Vtrap=[-0.0, -1194.0, -2047.0], Ctrap=1.0, frot=177.0, Vwall=0.10, wall_order=2)
        self.assertTrue(d.wz / (2 * pi) * .90 <= 1.865E6 <= d.wz / (2 * pi) * 1.10)

    def test_043015_axial_mode_consistency(self):
        """
        Checks to see if the axial frequency (which is identically equal to the COM mode eigenfrequency)
        agrees with the values obtained for varying frot in the lab on April 30th, 2015's data set.
        """
        a = ModeAnalysis(shells=10, Vtrap=[-0.0, -1750.0, -2000.0], Ctrap=1.0, frot=177.5, Vwall=2.13, wall_order=2)
        a.run()
        # Using units of hertz
        self.assertTrue(a.wz / (2 * pi) * .90 <= 1.576E6 <= 1.10 * a.wz / (2 * pi))
        # b = ModeAnalysis(shells=10 ,Vtrap=[-0.0,-1750.0,-2000.0], Ctrap = 1.0, frot=180.0, Vwall= 2.13, wall_order=2)
        # b.run_quiet()
        # self.assertTrue(b.wz/(2*pi)*.90<= 1.578E6 and 1.578E6 <= 1.10*b.wz/(2*pi))

    def test_051415_plane_stability(self):
        # Checks to make sure that we start to see a 1-2 plane transistion at appropriate frequencies.

        for w in np.linspace(187.5, 210, 20):
            a = ModeAnalysis(shells=10, Vtrap=[0.0,-1750.0,-2000.0], Ctrap = 1.0, frot=w, Vwall= 2.13, wall_order=2)
            a.run()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculationConsistency)
    unittest.TextTestRunner(verbosity=1).run(suite)

