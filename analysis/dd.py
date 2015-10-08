import scipy.constants as u
import numpy as np
import matplotlib.pylab as plt
import scipy.integrate
import scipy.interpolate


"""Functionality related to dynamical decoupling and noise."""


def filter_function_udd_spin_echo(w, tau, t_pi):
    """Spin Echo filter function
        Uys, H., PRL 103, 40501-40501 (2009)

    :param w: angular frequency  (rad/s)
    :param t: total decoupling duration (free-precession + pi-time)  (s)
    :param t_pi: pi-time  (s)
    :return: filter function value at t, w in range [0,1]
    """
    nz = 1/16  # normalization
    ftmp = (1 + np.exp(1j*w*tau) -
               2*np.exp(1j*w*tau/2.0)*np.cos(w*t_pi/2.0) )
    return nz*np.abs(ftmp*np.conj(ftmp))


def psd_from_vsd_func(w, vsd, eta):
    """ Compute interpolating function for power spectral density of
    of frequency fluctuations assuming a sensitivity of eta.

    :param w: angular frequency  (rad/sec)
    :param vsd: voltage spectral density  (V/sqrt(Hz))
    :param eta: pick-up coil sensitivity: b=-v/(eta*w)  (m**2)
    :return: interpf(w), an interpolating function for frequency power spectral
        density at w (rad/s) in units of (rad/s)**2/(rad/s)
    """
    uB = u.value("Bohr magneton")
    pi = np.pi
    g = 2.002
    hbar = u.hbar
    bsd = -1*vsd/(eta*w)  # magnetic field amplitude spectral density
    wsd = g*uB/(2*pi*hbar)*bsd  # frequency spectral density
    psd = wsd**2
    interpf = scipy.interpolate.interp1d(w, psd, kind=1)
    return interpf


def calc_coherence(tau, beta, ff, wmin, wmax):
    """ Compute spin coherence after free precession interval tau.
    Assume use of dynamical decoupling during tau.
        Uys, H., PRL 103, 40501-40501 (2009)

    :param tau: total free-precession time
            (including dynamical decoupling pulses)  (s)
    :param beta(w):  power spectral density of frequency
            fluctuations   (rad/s)**2/(rad/s)
    :param t_pi: pi-time for dynamical decoupling pulses  (s)
    :param ff(w): filter function
    :param wmin: minimum angular frequency (rad/s)
    :param wmax: maximum angular frequency (rad/s)
    :return: spin coherence [0, 1]
    """
    ws = np.linspace(wmin, wmax, 1000)
    integrand = [beta(w)*ff(w)/w**2 for w in ws]
    chi = 2.0/np.pi * scipy.integrate.simps(integrand, ws)
    return np.exp(-chi)