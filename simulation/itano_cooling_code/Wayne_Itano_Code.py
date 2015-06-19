__author__ = 'sbt'

# Began life as MATLAB code written by Charles Xu.
# Written for Python by Steven B. Torrisi, Summer 2015.

import time
import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.integrate as integ
import matplotlib.pyplot as plt
from scipy.constants import pi

from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

plt.rcParams['font.size'] = 16


class ItanoAnalysis:
    """
    Contains the tools to conduct calculations of total energy balance, torque, scattering events
    from a perpendicular doppler cooling beam acting incident on a Berylium ion plasma.
    """
    hbar = 1.054E-34
    c = 2.9979E8
    kB = 1.38E-23  # Boltzmann's constant
    lambd = 313.112E-9  # transition wavelength
    k = 2 * pi / lambd  # Be+ laser wavenumber; assuming offset dw is sufficiently
    m = 8.96 * 1.673E-27  # Be+ ion mass
    gamma0 = 2 * pi * 18.0E6  # natural linewidth
    r_bar = (hbar * k) ** 2 / (2 * m)  # recoil energy
    rnorm = r_bar / (hbar * k)  # A reduced recoil energy
    vk = gamma0 / (2 * k)  # Reduced line width
    pardiff=2*np.pi*9E6

    def __init__(self,
                 defaultoff=30.0E-6, defaultdet=-500E6, wr=2 * pi * 45.0E3, Tguess=1E-3, saturation=.5,
                 dens=2.77E9, ywidth=2.0E-6, radius=225.0E-6, quiet=True, spar=.2, wpar=2E6):

        """
        Define parameters which characterize the plasma and the incident doppler beam, such as the
        beam width in the y direction, radius and density, and saturation parameter for scattering interactions.

        :type defaultdet: float, default detunning frequency in Hz
        :type defaultoff: float, default doppler beam offset in meters
        :param wr: float, default angular frequency of rotation
        :param Tguess: float, default temperature guess
        :param saturation: float, saturation parameter for the beam-plasma interaction
        :param dens: float, density of plasma at (0,0)
        :param ywidth: float, describes the beam waist of the laser from a gaussian beam profile in the y-direction only
        :param radius: float, the radius of the plasma in consideration
        :param quiet: boolean, if false will print more output.
        """

        # Sets up parameters from initialization
        self.quiet = quiet  # controls verbosity of output; false for silent, true for more information
        self.d = defaultoff  # default offset of cooling laser
        self.det = defaultdet  # default detuning of cooling laser
        self.wr = wr  # angular rotation frequency in s^{-1}
        self.rp = radius  # radius of ion array in meters
        self.u0 = np.sqrt(2 * self.kB * Tguess / self.m)  # initial thermal velocity used in Boltzmann factor
        self.s0 = saturation  # saturation parameter for laser interaction (set to 0 to turn off saturation effects)
        self.sig0 = dens
        self.wy = ywidth


        # if you choose to include a parallel laser

        self.spar = spar
        self.wpar = wpar

        self.counter = 0

    def density(self, x, y):
        """
        Defines a density scaling which returns 0 outside of the circular bounds of the crystal.
        Useful because it allows for rectangular integration bounds (i.e. [-r,r])
        instead of complicated functional bounds which are difficult to work with.


        :rtype : float, either the height of the plasma in the z direction or 0 if out of bounds.
        :param x: x coordinate
        :param y: y coordinate
        """
        rad = x ** 2 + y ** 2
        if rad <= (self.rp ** 2):
            return np.sqrt(1 - (rad) / self.rp ** 2)
        else:
            return 0.0

    def dEavg(self, u, detun=None, offset=None):
        """

        Returns the average rate of change of energy from the perpendicular doppler cooling laser.

        :param u: float, which characterizes the Maxwell-Boltzmann distribution of velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb).
        :param detun: detuning frequency from the natural frequency of transistion. this should be negative.
        :param offset: offset in meters of the cooling laser from the center of the plasma.
        """
        if detun is None:
            detun = self.det
        if offset is None:
            offset = self.d

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr
        s0 = self.s0
        delta = 2. * detun / self.gamma0
        ret = integ.tplquad(lambda y, x, v:
                            self.density(x, y)
                            * np.exp(-(y - offset) ** 2 / wy ** 2) *
                            ((v) + 5 * self.rnorm / (6 * u)) * np.exp(-v ** 2) /
                            ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                              + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                            -np.inf, np.inf,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * rp, lambda x, y: rp)

    def dEavgAndParallel(self, u, detun=None, offset=None):
        """

        Returns the average rate of change of energy from the perpendicular doppler cooling laser.

        :param u: float, which characterizes the Maxwell-Boltzmann distribution of velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb).
        :param detun: detuning frequency from the natural frequency of transistion. this should be negative.
        :param offset: offset in meters of the cooling laser from the center of the plasma.
        """
        if detun is None:
            detun = self.det
        if offset is None:
            offset = self.d

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr
        s0 = self.s0
        delta = 2. * detun / self.gamma0

        alpha = (2 / 9) * self.hbar * self.k * rp ** 2 * np.pi ** (3 / 2) * self.spar  \
            /(s0*self.m*(1+2*self.spar+(2/self.gamma0)**2)*(self.pardiff)**2)

        ret = integ.tplquad(lambda y, x, v:
                            self.density(x, y)
                            * np.exp(-(y - offset) ** 2 / wy ** 2) *
                            ((v) + 5 * self.rnorm / (6 * u)) * np.exp(-v ** 2) /
                            ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                              + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                            -np.inf, np.inf,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * rp, lambda x, y: rp)

        # print("Integral Evaluated",ret[0],"with temperature",u**2*8.96*1.673E-27/(2*1.38E-23))
        self.counter += 1
        return u * ret[0] + alpha

    def totalscatter(self, ueq, detun=None, offset=None):
        """
        Returns the number of total scattering events per second from the perpendicular doppler cooling laser.


        :param ueq: float, which characterizes the Maxwell-Boltzmann distribution of velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb). Should be the eqilibrium temperature;
                run this after solving for a minimum U ideally.
        :param detun: detuning frequency from the natural frequency of transistion. this should be negative.
        :param offset: offset in meters of the cooling laser from the center of the plasma.
        """
        if detun is None:
            detun = self.det

        if offset is None:
            offset = self.d

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        delta = 2. * detun / self.gamma0
        constants = self.gamma0 * self.s0 * self.sig0 * ueq / np.sqrt(np.pi)
        ret = integ.tplquad(lambda y, x, v:
                            self.density(x, y)
                            * np.exp(-2 * (y - offset) ** 2 / wy ** 2) *
                            np.exp(-v ** 2) /
                            (1 + 2 * self.s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2) +
                             (delta - (wr * y / vk) - (ueq * v / vk)) ** 2),
                            -np.inf, np.inf,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * rp, lambda x, y: rp)

        return ret[0] * constants

    def totaltorque(self, ueq, detun=None, offset=None):
        """
            Returns the total torque imparted from the perpendicular doppler cooling laser.


        :param ueq: float, which characterizes the Maxwell-Boltzmann distribution of velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb). Should be the eqilibrium temperature;
                run this after solving for a minimum U ideally.
        :param detun: detuning frequency from the natural frequency of transistion. this should be negative.
        :param offset: offset in meters of the cooling laser from the center of the plasma.
        """

        if offset is None:
            offset = self.d

        if detun is None:
            detun = self.det

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        delta = 2. * detun / self.gamma0
        constants = (self.hbar * self.k) * self.gamma0 * self.s0 * self.sig0 * ueq / np.sqrt(np.pi)
        ret = integ.tplquad(lambda v, x, y:
                            self.density(x, y)
                            * np.exp(-2 * (y - offset) ** 2 / wy ** 2) * y *
                            np.exp(-v ** 2) /
                            (1 + 2 * self.s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2) +
                             (delta - (wr * y / vk) - (ueq * v / vk)) ** 2),
                            -1 * rp, rp,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * np.inf, lambda x, y: np.inf,
                            epsabs=1.49e-10, epsrel=1.49e-10)

        return ret[0] * constants

    """
    The following methods APPLY all of the methods above to facilitate numerical calculations.
    Currently supported are 'scans' along one dimension of detuning or offset,
    or a 2d scan of the detuning/offset parameter space.

    Additionally, some simple plot options will also be available (of June 11 2015)
    """

    def scan_detuning(self, offset=None):
        """
        Not ready yet.
        """
        if offset is None:
            offset = self.d
        return 0

    def scan_detuning_and_offset(self, detmin=-150.E6, detmax=-1.E6, detN=30,
                                 offmin=0, offmax=40.0E-6, offN=30,
                                 get_temp=True, get_torque=False, get_total_scatter=False, plot=False) -> list:
        """
            Conducts a two-dimensional scan across the detuning/offset parameter space.
            Returns an array of output elements. if plot is true, will generate plots for you of the temperature.


            :rtype : list which contains 3 elements:
                1) Teq element, a 2d array of the equilibrium temperatures, or the integer 0 if get_temp is false.
                2) Trq element, a 2d array of the net torques, or the integer 0 if get_torque is false.
                3) Sct element, a 2d array of the net scattering rate, or the integer 0 if get_total_scatter is false.

            :param detmin: Minimum detuning. Make negative, units of Hz.
            :param detmax: Maximum detuning. Make negative, units of Hz.
            :param detN:  Number of detunings to evaluate (interpolates linearly)
            :param offmin: Minimum offset. Make nonnegative, units of meters.
            :param offmax: Maximum offset. Make nonnegative, units of meters.
            :param offN:  Number of offsets to evaluate (interpolates linearly)
            :param get_temp: If true, will return temperature in the final array of results.
            :param get_torque: If true, will return the torques in the final array of results.
            :param get_total_scatter: If true, will return the scattering rates in the final array of results.
            :param plot:
            """
        Tmax = 1.0E-1
        Tmin = 1.0E-5
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        Trq = 0
        Sct = 0
        df = np.linspace(detmin, detmax, detN)
        doff = np.linspace(offmin, offmax, offN)
        dw = 2 * pi * df
        U = np.zeros([len(df), len(doff)])
        if get_torque is True:
            Trq = np.zeros([len(df), len(doff)])
        if get_total_scatter is True:
            Sct = np.zeros([len(df), len(doff)])
        for W in range(len(dw)):
            for D in range(len(doff)):
                if self.quiet is False:
                    print("-------------")
                    print("Feeding in detuning/offset", dw[W] / (2 * 3.14159), doff[D], W * (len(doff)) + D + 1, "of",
                          len(df) * len(doff), "Evaluations",
                          100 * (W * (len(doff)) + D + 1) / (len(df) * len(doff)), "%")

                # The difference between each of these three tries is that
                try:
                    U[D][W] = opt.brentq(self.dEavg, umin, umax, args=(dw[W], doff[D]), xtol=1e-4, rtol=3.0e-7)
                except:
                    try:
                        U[D][W] = opt.brentq(self.dEavg, umin * .1, umax * 10, args=(dw[W], doff[D]), xtol=1e-4,
                                             rtol=3.0e-7)
                    except:

                        try:
                            U[D][W] = opt.brentq(self.dEavg, umin * .01, umax * 100, args=(dw[W], doff[D]), xtol=1e-4,
                                                 rtol=3.0e-7)
                        except:
                            U[D][W] = 313.0  # Obviously an error occured

                if get_torque:
                    Trq[D][W] = self.totaltorque(U[D][W], dw[W], doff[D])

                if get_total_scatter:
                    Sct[D][W] = self.totalscatter(U[D][W], dw[W], doff[D])

                if self.quiet is False:
                    print("u/T recorded:", U[D][W], ((U[D][W]) * U[D][W] * 8.96 * 1.673E-27 / (2 * 1.38E-23)), "after",
                          self.counter, "integral calls.")

                self.counter = 0
        if get_temp is True:
            Teq = np.array([u ** 2 * self.m / (2 * self.kB) for u in U])
        else:
            Teq = 0
        return [Teq, Trq, Sct]

    def getueq(self, detun=None, offset=None):
        if detun is None:
            detun = self.det
        if offset is None:
            offset = self.d

        Tmax = 1.0E-1
        Tmin = 1.0E-5
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        ret = opt.brentq(self.dEavg, umin, umax, args=(detun, offset), xtol=1e-4, rtol=3.0e-7)
        return ret[0]

    """
    def getrootdetun(self,detun):
        global GAMMA, VK, RNORM, U0, RP, WR, WY, d,S0,SIG
        print("Trying detuning",detun)
        return (opt.brentq(dEver3,umin,umax,args=(detun,3.98794261e-05),xtol=1e-6, rtol=3.0e-8))**2*m/(2*kB)
    """


if __name__ == '__main__':
    a = ItanoAnalysis(quiet=False)

    A = a.scan_detuning_and_offset(detN=3, offN=3, get_torque=True, get_total_scatter=True)
    df = np.linspace(-150.0E6, -1E6, 3)
    doff = np.linspace(0, 40.0E-6, 3)
    print(A[0], A[1], A[2])
    plt.rcParams['font.size'] = 16
    plt.figure()
    CS = plt.pcolor(df * 1.0E-6, 1.0E6 * doff, A[0] * 1000)
    plt.show()
    plt.figure()
    CS2 = plt.pcolor(df * 1.0E-6, 1.0E6 * doff, A[1])
    plt.show()
    """
    Tmax = 1.0E-1
    Tmin = 1.0E-5
    umax = np.sqrt(Tmax * kB * 2 / m)
    umin = np.sqrt(Tmin * kB * 2 / m)
    global counter
    counter = 0
    print("Beginning...")

    # opt.fmin(getroot,[-50.0E6,15.0E-6],args=(detun))
    # print(opt.minimize(getroot,[-50.0E6,50.0E-6],method='TNC',bounds=[[-160.0E6,-1.0E6],[0,220.0E-6]]))

    # print(opt.minimize(getrootdetun,[-10.0E6] ,method='TNC',bounds=[[-150.0E6,-5.0E6]]))
    # return


    # plt.figure(figsize=(12, 10))
    # b=plt.plot(df,Teq)
    # b=plt.xlabel('laser freq. detuning (Hz)')
    # b=plt.ylabel('equilibrium temp (K)')
    # b=plt.title('Plot of the Equilibrium Temperature')
    plt.figure(figsize=(12, 10))
    # cmap = plt.get_cmap('PiYG')
    # levels = MaxNLocator(nbins=15).tick_values(Teq.min(), Teq.max())
    # norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    # im = plt.pcolormesh(df, doff, Teq, cmap=cmap, norm=norm)
    CS = plt.contour(doff, df, Teq)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('Simplest default with labels')
    a=plt.plot(df,Teq)
    a=plt.xlabel('laser freq. detuning (Hz)')
    a=plt.ylabel('equilibrium temp (K)')
    a=plt.title('Plot of the Equilibrium Temperature')
    a=plt.yscale('log', nonposy='clip')

    returnplt.figure(figsize=(12, 10))
    c=plt.plot(df,Sc)
    c=plt.xlabel('laser freq. detuning (Hz)')
    c=plt.ylabel('T`otal Scattering Rate')
    c=plt.title('Detuning vs Total Scattering Rate')

    plt.figure(figsize=(12, 10))
    d=plt.plot(df,To)
    d=plt.xlabel('laser freq. detuning (Hz)')
    d=plt.ylabel('Net Torque by Laser')
    d=plt.title('Detuning vs Net Torque by Laser')

    print(U)
    """
