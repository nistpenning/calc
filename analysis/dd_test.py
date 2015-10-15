import unittest
import importlib

import scipy.constants as u
import numpy as np
import matplotlib.pylab as plt
import scipy.integrate

from dd import *

debug = True


class TestFilterFunction(unittest.TestCase):

    def test_se_ff(self):
        # expectation: if t_arm=1ms, peak in sensitivity is at 500 kHz;
        # analytic form holds for t_pi = 0
        arm_t = 0.5e-3
        t_pi = 0
        fd = np.linspace(0, 5000, 100)
        filter = filter_function_udd_spin_echo(w=2*np.pi*fd, tau=2*arm_t+t_pi,
                                               t_pi=t_pi)
        analytic = 16*(np.sin(2*np.pi*fd*2*arm_t/4))**4
        if debug:
            plt.clf()
            plt.plot(fd, filter, 'y-', alpha=0.5)
            plt.plot(fd, analytic, '--g', alpha=0.5)
            plt.title("ff for t_arm = {}".format(arm_t))
            plt.ylabel("filter pass")
            plt.xlabel("Freq [Hz]")
            plt.show()
        residual = analytic-filter
        self.assertTrue(sum(residual) < 1e-3)


class TestCalculatedCoherence(unittest.TestCase):
    """
    See dd_test.lyx for the approach used here.
    """
    def constant_vsd(self, ws, v0):
        """Constant voltage spectral density

        :param v0: amplitude  (v)
        :return: voltage spectral density  (v/sqrt(rad/s))
        """
        return np.array([v0 for w in ws])

    def coh_for_constant_vsd(self, t, eta, v0):
        """Expected coherence for VSD test_constant_vsd(). This is a
        consistency check.
        """
        pi = np.pi
        hbar = u.hbar
        uB = u.value("Bohr magneton")
        c1 = (2/pi)*(2*uB/(2*pi*u.hbar*eta))**2*(pi*16/192)
        return np.exp(-1*c1*v0**2*t**3)

    def test_coh(self):
        pi = np.pi
        taus = np.linspace(1e-3, 0.1, 20)
        wmin = 0.1*2*pi
        wmax = 1e3*2*pi
        v0 = 5e-8
        eta = 1
        ws = np.linspace(wmin, wmax, 3000)
        vsd = self.constant_vsd(ws=ws, v0=v0)
        beta = psd_from_vsd_func(ws, vsd, eta=eta)
        coh = np.zeros(len(taus))
        for i, tau in enumerate(taus):
            ff_at_tau = lambda w: filter_function_udd_spin_echo(w, tau, t_pi=0)
            coh[i] = calc_coherence(tau, beta=beta, ff=ff_at_tau,
                                wmin=wmin, wmax=wmax)
        coh_expected = [self.coh_for_constant_vsd(tau, eta=eta, v0=v0) for tau in taus]
        mean_residual = np.mean(coh_expected - coh)
        if debug:
            plt.plot(taus, coh, 'k+')
            plt.plot(taus, coh_expected, 'r-')
            plt.ylim(-0.1, 1.1)
            plt.show()
        self.assertTrue(mean_residual < 1e-6)

if __name__ == '__main__':
    unittest.main(exit=False)












