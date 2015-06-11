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

plt.rcParams['font.size'] = 16


class ItanoAnalysis:
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

    def __init__(self, detmin=-150.E6, detmax=-1.E6, detN=30, offmin=0, offmax=40.0E-6, offN=30,
                 defaultoff=0, defaultdet=-20E-6, wr=2 * pi * 45.0E3, Tguess=1E-3, saturation=.5,
                 dens=2.77E9, ywidth=2.0E-6, quiet=True):
        self.quiet = quiet
        self.df = np.linspace(detmin, detmax, detN)  # Interpolate linearly from lowest df to highest in N steps
        self.dw = 2 * pi * self.df
        self.d = defaultoff
        self.det = defaultdet
        self.doff = np.linspace(offmin, offmax, offN)
        d = 40.00E-6  # offset of laser beam from the center of the array in meters,
        # positive offsets corresponds to miving the laser beam to the side of
        # the cloud (or array) that is receding from the laser beam direction
        self.wr = wr  # angular rotation frequency in s^{-1}
        self.rp = 225.0E-6  # radius of ion array in meters
        self.u0 = np.sqrt(2 * self.kB * Tguess / self.m)  # initial thermal velocity used in Boltzmann factor
        self.s0 = saturation
        self.sig0 = dens
        self.wy = ywidth

        self.counter = 0

    def density(self, x, y):
        """
        Defines a density scaling which returns 0 outside of the circular bounds of the crystal.
        Useful because it allows for rectangular integration bounds (i.e. [-r,r])
        instead of complicated functional bounds which are difficult to work with.
        """
        rad = x ** 2 + y ** 2
        if rad <= (self.rp ** 2):
            return np.sqrt(1 - (rad) / self.rp ** 2)
        else:
            return 0

    def dEnosaturation(self, u, detun, offset=None):

        if offset is None:
            offset = self.d

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr
        delta = 2. * detun / self.gamma0
        ret = integ.tplquad(lambda y, x, v: \
                                self.density(x, y) \
                                * np.exp(-(y - offset) ** 2 / wy ** 2) * \
                                ((v) + 2 * self.rnorm / u) * np.exp(-v ** 2) / \
                                (1 + (delta - (wr * y / vk) - (u * v / vk)) ** 2), \
                            -np.inf, np.inf, \
                            lambda x: -1 * rp, lambda x: rp, \
                            lambda x, y: -1 * rp, lambda x, y: rp)

        # print("Integral Evaluated",ret[0],"with temperature",u**2*8.96*1.673E-27/(2*1.38E-23))
        self.counter += 1
        return ret[0]

    def dEavg(self, u, detun, offset=None):

        if offset is None:
            offset = self.d

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        delta = 2. * detun / self.gamma0
        ret = integ.tplquad(lambda y, x, v: \
                                self.density(x, y) \
                                * np.exp(-(y - offset) ** 2 / wy ** 2) * \
                                ((v) + 2 * self.rnorm / u) * np.exp(-v ** 2) / \
                                ((1 + 2 * S0 * np.exp(-2 * (y - d) ** 2 / wy ** 2) \
                                  + (delta - (wr * y / vk) - (u * v / vk)) ** 2)), \
                            -np.inf, np.inf, \
                            lambda x: -1 * rp, lambda x: rp, \
                            lambda x, y: -1 * rp, lambda x, y: rp)

        # print("Integral Evaluated",ret[0],"with temperature",u**2*8.96*1.673E-27/(2*1.38E-23))
        self.counter += 1
        return ret[0]

    def totalscatter(self, ueq, detun, offset=None):

        if offset is None:
            offset = self.d
        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        delta = 2. * detun / self.gamma0
        constants = self.gamma0 * self.s0 * self.sig0 * ueq / np.sqrt(np.pi)
        ret = integ.tplquad(lambda y, x, v: \
                                \
                                self.density(x, y) \
                                * np.exp(-2 * (y - offset) ** 2 / wy ** 2) * \
                                np.exp(-v ** 2) / \
                                (1 + 2 * S0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2) + \
                                (delta - (wr * y / vk) - (ueq * v / vk)) ** 2), \
                                -np.inf, np.inf, \
                            lambda x: -1 * rp, lambda x: rp, \
                            lambda x, y: -1 * rp, lambda x, y: rp)

        return ret[0] * constants

    def totaltorque(self, ueq, detun = None, offset=None):

        if offset is None:
            offset = self.d

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        delta = 2. * detun / self.gamma0
        constants = (self.hbar * self.k) * self.gamma0 * self.s0 * self.sig0 * ueq / np.sqrt(np.pi)
        ret = integ.tplquad(lambda v, x, y: \
                            self.density(x, y) \
                            * np.exp(-2 * (y - offset) ** 2 / wy ** 2) * y * \
                            np.exp(-v ** 2) / \
                            (1 + 2 * self.s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2) + \
                             (delta - (wr * y / vk) - (ueq * v / vk)) ** 2), \
                            -1 * rp, rp, \
                            lambda x: -1 * rp, lambda x: rp, \
                            lambda x, y: -1 * np.inf, lambda x, y: np.inf, \
                            epsabs=1.49e-10, epsrel=1.49e-10)

        return ret[0] * constants

    def scan_detuning(self, offset=None):

        return 0

    def scan_detuning_and_offset(self, get_temp = True, get_torque=False, get_total_scatter=False):

        Tmax = 1.0E-1
        Tmin = 1.0E-5
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        Trq = 0
        Sct = 0
        df = np.linspace(-80.0E6, -5.0E6, 40)
        doff = np.linspace(0.00, 40.0E-6, 40)
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
                        len(df) * len(doff), "Evaluations", \
                        100 * (W * (len(doff)) + D + 1) / (len(df) * len(doff)), "%")

                #The difference between each of these three tries is that
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

                if self.quiet is False:
                    print("u/T recorded:", U[D][W], ((U[D][W]) * U[D][W] * 8.96 * 1.673E-27 / (2 * 1.38E-23)), "after",
                        self.counter, "integral calls.")

                self.counter = 0

        if get_temp is True:
                Teq = np.array([u ** 2 * self.m / (2 * self.kB) for u in U])
        else:
            Teq = 0


        return [Teq, Trq, Sct]




    """
    def getroot(self,detoffset):
        global GAMMA, VK, RNORM, U0, RP, WR, WY, d,S0,SIG
        print("Trying detuning",detoffset[0]," and offset ",detoffset[1])
        return (opt.brentq(dEver3,umin,umax,args=(detoffset[0],detoffset[1]),xtol=1e-6, rtol=3.0e-8))**2*m/(2*kB)

    def getrootdetun(self,detun):
        global GAMMA, VK, RNORM, U0, RP, WR, WY, d,S0,SIG
        print("Trying detuning",detun)
        return (opt.brentq(dEver3,umin,umax,args=(detun,3.98794261e-05),xtol=1e-6, rtol=3.0e-8))**2*m/(2*kB)
    """

if __name__ == '__main__':
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
