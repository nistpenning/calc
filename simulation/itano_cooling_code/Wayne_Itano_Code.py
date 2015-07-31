__author__ = 'sbt'

""" Began life as MATLAB code written by Charles Xu.
    Written for Python by Steven B. Torrisi, Summer 2015.
    Contains a class called ItanoAnalysis, based primarily on the PRA Paper by W. Itano from 1988
    "Perpendicular laser cooling of a rotating ion plasma in a Penning trap"
    with modifications which take into account a rotating wall potential, as well as
    some recoil heating from a parallel laser beam frmo the PRA paper by W. Itano from 1982
    "Laser Cooling of Atoms Stored in Harmonic and Penning Traps".

    July 21 Revision

"""
import numpy as np
import scipy.optimize as opt
import scipy.integrate as integ
import matplotlib.pyplot as plt
from scipy.constants import pi

plt.rcParams['font.size'] = 20


class ItanoAnalysis:
    """
    Contains the tools to conduct calculations of total energy balance, torque, scattering events
    from a perpendicular doppler cooling beam acting incident on a Berylium ion plasma.

    Methods inside the function allow one to find, for a plasma of a given radius and density,
    the thermal equilibrium which results from a perpendicular doppler cooling laser acting on the
    plasma.
    """
    # Define fundamental constants
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

    def __init__(self,
                 defaultoff=30.0E-6, defaultdet=-50E6, Tguess=1E-3,
                 wr=2 * pi * 45.0E3,
                 saturation=.5, dens=2.77E9, ywidth=2.0E-6, radius=225.0E-6, quiet=True,
                 parheating=False, spar=.2, wpar=-gamma0/2):

        """
        Define parameters which characterize the plasma and the incident doppler beam, such as the
        beam width in the y direction, radius and density, and saturation parameter for scattering interactions.

        :type defaultdet: float, default detunning frequency in Hz
        :type defaultoff: float, default doppler beam offset in meters
        :param wr: float, default angular frequency of rotation
        :param Tguess: float, default temperature guess in Kelvin
        :param saturation: float, dim.less saturation parameter for the beam-plasma interaction
        :param dens: float, density of plasma at (0,0)
        :param ywidth: float, describes the beam waist of the laser from a gaussian beam profile
                                in the y-direction only
        :param radius: float, the radius of the plasma in consideration
        :param quiet: boolean, if false will print more output.

        :param parheating: boolean, if on, will add parallel laser heating to dE calculations
        :param spar: float, Saturation parameter for parallel cooling laser
        :param wpar: float, Detuning for the parallel cooling laser
        """

        # Store parameters from initialization
        self.quiet = quiet  # controls verbosity of output; false for silent,
        #                       true for more information
        self.d = defaultoff  # default offset of cooling laser
        self.det = defaultdet  # default detuning of cooling laser
        self.wr = wr  # angular rotation frequency in s^{-1}
        self.rp = radius  # radius of ion array in meters
        self.u0 = np.sqrt(
            2 * self.kB * Tguess / self.m)  # initial thermal velocity used in Boltzmann factor
        self.s0 = saturation  # saturation parameter for laser interaction (set to 0 to
        #                                turn off saturation effects)
        self.sig0 = dens
        self.wy = ywidth

        # if you choose to include a parallel laser

        self.parheating = parheating

        self.spar = spar
        self.wpar = wpar
        self.alpha = (1 / 9) * self.hbar * self.k * self.rp ** 2 * np.pi ** (3 / 2) * self.spar \
                     / (self.s0 * self.m * (1 + self.spar))

        self.counter = 0

    def density(self, x, y, radius=None):
        """
        Defines a density scaling which returns 0 outside of the circular bounds of the plasma.
        Useful because it allows for rectangular integration bounds (i.e. [-r,r])
        instead of complicated functional bounds which are difficult to work with
        (such as [-\sqrt{1-x^2}, \sqrt{1-x*2}] for example).

        Uses the plasma radius as defined in the class (self.rp)

        :rtype : float: either the height of the plasma in the z direction or 0 if out of bounds.
        :param x: x coordinate
        :param y: y coordinate
        """
        rad = x ** 2 + y ** 2

        if radius is None:
            radius = self.rp
        if rad <= (radius ** 2):
            return np.sqrt(1 - rad / radius ** 2)
        else:
            return 0.0

    def dEavg(self, u, detun=None, offset=None):
        """

        Returns the average rate of change of energy from the perpendicular doppler cooling laser.

        :param u: float, which characterizes the Maxwell-Boltzmann distribution of velocity for
                    the ions.
                Temperature can be inferred by u^2 *m /(2kb).
        :param detun: detuning frequency from the natural frequency of transistion.
                        This should be negative.
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
                            * np.exp(-2 * (y - offset) ** 2 / wy ** 2) *
                            (v + 5 * self.rnorm / (3 * u)) * np.exp(-v ** 2) /
                            ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                              + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                            -np.inf, np.inf,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * rp, lambda x, y: rp)
        self.counter += 1
        if self.parheating is True:
            return ret[0] * u + self.alpha
        if self.parheating is False:
            return ret[0]

    def dEavg_fine(self, u, detun=None, offset=None,num=50):
        """

            Returns the average rate of change of energy from the perpendicular
            doppler cooling laser.

            :param u: float, which characterizes the Maxwell-Boltzmann distribution of
                    velocity for the ions.
                    Temperature can be inferred by u^2 *m /(2kb).
            :param detun: detuning frequency from the natural frequency of transistion.
                    This should be negative.
            :param offset: offset in meters of the cooling laser from the center of the plasma.
            :param num: integeer which describes the precision of the calculation
                        (number of subdivisions)
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
        # Ymin and Ymax chosen for intensities of < one thousandth of the peak
        ymin = offset - self.wy * 2
        ymax = offset + self.wy * 2
        yrange = np.linspace(ymin, ymax, num)
        ystep = yrange[1] - yrange[0]

        dE = 0

        # Integrate over small pieces of the y-range of the crystal.
        for Y in yrange:
            ret = integ.tplquad(lambda y, x, v:
                                self.density(x, y)
                                * np.exp(-2 * (y - offset) ** 2 / wy ** 2) *
                                ((v) + 5 * self.rnorm / (3 * u)) * np.exp(-v ** 2) /
                                ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                                  + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                                -np.inf, np.inf,
                                lambda x: -1 * rp, lambda x: rp,
                                lambda x, y: Y, lambda x, y: Y + ystep)
            # Notice the range of the y integrations change in minute ways
            # from call to call.
            print("From range [", Y, Y + ystep, "]:", ret[0])
            dE += ret[0]
        self.counter += 1
        if self.parheating is True:
            return dE * u + self.alpha
        if self.parheating is False:
            return dE

    def dEavg_no_rotating_wall(self, u, detun=None, offset=None):
        """
        NOT SURE IF THIS IS CURRENTLY WORKING! As of July 16
        Returns the average rate of change of energy from the
        perpendicular doppler cooling laser, without the rotating wall.

        Uses the same methodology as the fine integration.

        :param u: float, which characterizes the Maxwell-Boltzmann distribution of velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb).
        :param detun: detuning frequency from the natural frequency of transistion. this should be negative.
        :param offset: offset in meters of the cooling laser from the center of the plasma.
        """
        # print("loud and clear")
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


        # print('The terms in question:')
        # print(u*mult)
        # print(mult*wr*1.0E-6/(self.k*u))
        # print(mult*5 *self.hbar*self.rnorm / (3*u))

        # print(mult*wr/self.k)
        # print(5*self.rnorm*self.hbar /3)

        # print("denominator terms:")
        # print("wr*y/vk ~~",wr*1E-6/vk)
        # print("u*v/vk~~",u/vk)
        """
        ret = integ.tplquad(lambda x, y, v:
                        mult*self.density(x, y)
                        * np.exp(-2*(y - offset) ** 2 / wy ** 2) *
                        (v*self.hbar + wr*y/(self.k*u) + 5*self.rnorm*self.hbar / (3*u))
                        * np.exp(-v ** 2) /
                        ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                        + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                        -np.inf, np.inf,
                        lambda x: -1 * rp, lambda x: rp,
                        lambda x, y: -1 * rp, lambda x, y: rp)
        """
        print("Evaluating")
        ymin = offset - self.wy * 1.9
        ymax = offset + self.wy * 1.9
        mult = 1.0

        num = 50
        yrange = np.linspace(ymin, ymax, num)
        ystep = yrange[1] - yrange[0]
        dE = 0
        for Y in yrange:
            ret = integ.tplquad(lambda y, x, v:
                                self.density(x, y)
                                * np.exp(-2 * (y - offset) ** 2 / wy ** 2) *
                                (mult * v * u * self.hbar + mult * wr * y / (
                                    self.k) + mult * 5 * self.rnorm * self.hbar / 3)
                                * np.exp(-v ** 2) /
                                ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                                  + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                                -np.inf, np.inf,
                                lambda x: -1 * rp, lambda x: rp,
                                lambda x, y: Y, lambda x, y: Y + ystep)
            dE += ret[0]
        self.counter += 1
        if self.parheating is True:
            return dE * u / mult + self.alpha
        if self.parheating is False:
            return dE / mult

    def totalscatterperp(self, ueq, detun=None, offset=None):
        """
        Returns the number of total scattering events per second from the perpendicular
                doppler cooling laser.

        Warning; This evaluates the scattering rate over one call of the triple integral
        function and may be subject to errors under certain conditions. Please consider using
        the peaktopeaktorque method, which computes the scattering in a fine way.

        :param ueq: float, which characterizes the Maxwell-Boltzmann distribution of
                velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb). Should be the eqilibrium temperature;
                run this after solving for a minimum U ideally.
        :param detun: detuning frequency from the natural frequency of transistion.
                        This should be negative.
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
        constants = self.gamma0 * self.s0 * self.sig0 / np.sqrt(np.pi)
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

        WARNING; This evaluates the torque rate over one call of the triple integral
        function and has been clearly demonstrated to have errors under certain
        conditions. Please consider using
        the peaktopeaktorque method, which computes the torque in a fine way.

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
        constants = (self.hbar * self.k) * self.gamma0 * self.s0 * self.sig0 / np.sqrt(np.pi)
        ret = integ.tplquad(lambda v, x, y:
                            self.density(x, y)
                            * np.exp(-2 * (y - offset) ** 2 / wy ** 2) * y *
                            np.exp(-v ** 2) /
                            (1 + 2 * self.s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2) +
                             (delta - (wr * y / vk) - (ueq * v / vk)) ** 2),
                            -1 * rp, rp,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * np.inf, lambda x, y: np.inf,
                            epsabs=1.49e-15, epsrel=1.49e-15)

        return ret[0] * constants

    def scattering_and_torque_fine(self, ueq, detun=None, offset=None):
        """
            Returns the maximum difference in torque (at differing points of y)
            imparted from the perpendicular doppler
            cooling laser.

            Returns the scattering, torque, and peak-to-peak torque intensity.


        :param ueq: float, which characterizes the Maxwell-Boltzmann distribution of velocity for the ions.
                Temperature can be inferred by u^2 *m /(2kb). Should be the eqilibrium temperature;
                run this after solving for a minimum U ideally.
        :param detun: detuning frequency from the natural frequency of transistion. this should be negative.
        :param offset: offset in meters of the cooling laser from the center of the plasma.
        """

        if np.isnan(ueq):
            return [np.nan, np.nan, np.nan]

        if offset is None:
            offset = self.d

        if detun is None:
            detun = self.det

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        delta = 2. * detun / self.gamma0

        constants = (self.hbar * self.k) * self.gamma0 * self.s0 * self.sig0 / np.sqrt(np.pi)

        # Only interested in Intensities above one thousandth of the peak:
        # 1.85 comes from
        # log .001 = -6.9
        # so e^-6.9 = .001
        # to solve for the offset at which intensity is at one thousandth,
        # -2 (y-d)^2 /wy^2 = -6.9
        # means (y-d)^2 = 6.9/2 * wy^2
        # y-d = 1.85 * wy
        # So we concern ourselves with
        # |y| <= 1.85*wy +d
        # Round up to 1.9 for good measure

        maxy = offset + self.wy * 1.9
        miny = offset - self.wy * 1.9
        # maxy = rp
        # miny = -rp

        num = 150

        yrange = np.linspace(miny, maxy, num=num)
        ystep = yrange[1] - yrange[0]

        torques = np.zeros(num)
        scattering = np.zeros(num)
        for i in range(len(yrange)):

            if ueq != np.NAN:
                y = yrange[i]
                xmin = -np.sqrt(rp - y ** 2)
                xmax = np.sqrt(rp - y ** 2)

                scattering[i] = integ.dblquad(lambda v, x:
                                              self.density(x, y)
                                              * np.exp(-2 * (y - offset) ** 2 / wy ** 2) *
                                              np.exp(-v ** 2) /
                                              (1 + 2 * self.s0 * np.exp(
                                                  -2 * (y - offset) ** 2 / wy ** 2) +
                                               (delta - (wr * y / vk) - (ueq * v / vk)) ** 2),
                                              xmin, xmax,
                                              lambda v: -1 * np.inf, lambda v: np.inf,
                                              epsabs=1.49e-7, epsrel=1.49e-7)[0]
            else:
                scattering[i] = 0

        totalscatter = constants / (self.k * self.hbar) * np.sum(scattering) * ystep
        # print("By the way the net torque is",np.sum(torques)*ystep)

        nettorque = integ.trapz(y=constants * scattering * yrange, dx=ystep)

        peak2peak = np.max(torques) - np.min(torques)
        return [totalscatter, nettorque, peak2peak]

    def scan_detuning(self, offset=None):
        """
        Not ready yet. Will allow you to generate a one-dimensional scan of detuning.
        """
        if offset is None:
            offset = self.d
        return 0

    def scan_detuning_and_offset(self, detmin=-150.E6, detmax=-1.E6, detN=30,
                                 offmin=0, offmax=40.0E-6, offN=30,
                                 get_scatter_and_torque=False,
                                 plot=False, no_rotating_wall=False,
                                 detuninglist=None, offsetlist=None,
                                 temphighprecision=False) -> list:
        """
            Conducts a two-dimensional scan across the detuning/offset
            parameter space.
            Returns an array of output elements.

            :rtype : list which contains 3 elements:
                1) Teq element, a 2d array of the equilibrium temperatures.
                2) Trq element, a 2d array of the net torques, or the integer
                    0 if get_scatter_and_torque is false.
                3) Sct element, a 2d array of the net scattering rate, or the
                    integer 0 if get_scatter_and_torque is false.
                4) P2p element, a 2d array of the peak-to-peak torque

            :param detmin: Minimum detuning. Make negative, units of Hz.
            :param detmax: Maximum detuning. Make negative, units of Hz.
            :param detN:  Number of detunings to evaluate
                            (interpolates linearly)
            :param offmin: Minimum offset. Make nonnegative, units of meters.
            :param offmax: Maximum offset. Make nonnegative, units of meters.
            :param offN:  Number of offsets to evaluate (interpolates linearly)
            :param get_torque: If true, returns the torques in the result array.
            :param get_total_scatter: If true, returns the scattering rates.
            :param detuninglist: Allows for the use of a custom list of
                            detunings instead of a linear interpolation.
            :param offsetlist: Allows for the use of a custom list of offsets instead
                                        of a linear interpolation.
            :param parheating: boolean, turns on/off parallel heating.
            :param highprecision: boolean, turns the equilibrium temperature calculation
                                to high precision mode.
            :param plot: Not implemented, but could be used to generate plots automatically.
            """
        Tmax = 1.0E-1
        Tmin = 1.0E-5
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        Trq = 0
        Sct = 0
        P2p = 0  # Peak to peak torque

        # Goes down a tree to decide which dEfunction to use.
        # Notice that certain functions overrule other ones. For example,
        # dEavg_fine overrules DEavgAndParallel. One must simply implement
        # new functions if they wish to have certain functions called with a variety of
        # features.

        dEfunction = self.dEavg

        if no_rotating_wall is True:
            dEfunction = self.dEavg_no_rotating_wall

        if temphighprecision is True:
            dEfunction = self.dEavg_fine

        if detuninglist is None:
            df = np.linspace(detmin, detmax, detN)
        else:
            df = detuninglist

        if offsetlist is None:
            doff = np.linspace(offmin, offmax, offN)
        else:
            doff = offsetlist

        dw = 2 * pi * df

        U = np.zeros([len(doff), len(df)])

        if get_scatter_and_torque is True:
            Trq = np.zeros([len(doff), len(df)])
            Sct = np.zeros([len(doff), len(df)])
            P2p = np.zeros([len(doff), len(df)])

        for W in range(len(df)):
            for D in range(len(doff)):
                if self.quiet is False:
                    print("-------------")
                    print("Feeding in detuning/offset", dw[W] / (2 * 3.14159),
                          doff[D], W * (len(doff)) + D + 1, "of",
                          len(df) * len(doff), "Evaluations. %2.2f" % (100 *
                                                                       (W * (len(doff)) + D + 1) / (
                                                                           len(df) * len(doff))),
                          "%")

                try:
                    U[D][W] = opt.brentq(dEfunction, umin, umax,
                                         args=(dw[W], doff[D]), xtol=1e-4, rtol=3.0e-7)
                except:
                    Nlist = [10 ** 2, 10 ** 4, 10 ** 6, 10 ** 8]
                    succeed = False
                    for factor in Nlist:
                        if succeed:
                            break
                        try:
                            # print("upping max...")
                            # print("Trying with umax",uMAX)
                            U[D][W] = opt.brentq(dEfunction, umin / 10, umax * factor,
                                                 args=(dw[W], doff[D]), xtol=1e-5, rtol=3.0e-8)
                            succeed = True
                        except:
                            pass

                    if not succeed:
                        U[D][W] = np.nan
                        print('Warning! The integral could not evaluate for some reason, '
                              'and was stored as NaN.')

                if get_scatter_and_torque:
                    results = self.scattering_and_torque_fine(U[D][W], dw[W], doff[D])
                    Sct[D][W] = results[0]
                    Trq[D][W] = results[1]
                    P2p[D][W] = results[2]

                if self.quiet is False:
                    print("u/T recorded:", U[D][W],
                          ((U[D][W]) * U[D][W] * 8.96 * 1.673E-27 / (2 * 1.38E-23)), "after",
                          self.counter, "integral calls.")

                self.counter = 0

        Teq = np.array([u ** 2 * self.m / (2 * self.kB) for u in U])

        return [Teq, Trq, Sct, P2p]

    def getueq(self, dEfunction=None, detun=None, offset=None, returntemperature=False):
        """
        Returns the equilibrium Maxwell-Boltzmann parameter for the planar velocity.
        or, if you set returntemperature to True, it will return the temperature isntead.

        :type detun: float: The detuning to use for the perpendicular cooling laser.
        :type offset: float: The offset from the center of the plasma to use for the perpendicular cooling laser.
        """
        if detun is None:
            detun = self.det
        if offset is None:
            offset = self.d
        if dEfunction is None:
            dEfunction = self.dEavg

        Tmax = 1.0E-1
        Tmin = .40E-3
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        try:
            ret = opt.brentq(dEfunction, umin, umax, args=(detun, offset), xtol=1e-4,
                             rtol=3.0e-7)

        except:
            for power in [10,1000,1E5,1E7]:
                try:
                    ret = opt.brentq(dEfunction, umin, umax*power,
                                     args=(detun, offset), xtol=1e-4,rtol=3.0e-7)
                    break
                except:
                    ret=np.nan


        if returntemperature is False:
            return ret[0]
        else:
            return ret[0] ** 2 * self.m / (2 * self.kB)

    def plot_temperature_result(self, temp, df, doff, showminimum=True,
                                title="Equilibrium Temperature (mK)", contourlevels=None,
                                maxtemp=None, mintemp=None):
        """
        Input the 0th array which returns from scan detuning and offset to
        make a plot of the temperature.

        :param temp: Array of temperatures to plot
        :param df: x-axis: detuning in Megahertz
        :param doff: y-axis: offset in micron
        :param showminimum: boolean, if true plots a star at the minimum temperature
        :param title: string, titles the plot
        :param contourlevels: list of values to plot temperature contours at
        :param maxtemp: maximum temperature to cut off temperatures above, and to put
            the top of the color plot at
        :param mintemp: minimum temperature to ground the color range at
        Assumptions about input parameters are as follows:


        """

        plt.rcParams['font.size'] = 22
        Teq = temp * 1000

        Teqmask = np.ma.array(Teq, mask=np.isnan(Teq))
        if maxtemp is not None:
            Teqmask = np.ma.array(Teqmask, mask=Teq >= maxtemp)

        if maxtemp is None:
            maxtemp = np.nanmax(Teqmask)
        if mintemp is None:
            mintemp = np.nanmin(Teqmask)
        offsetmin, freqmin = np.unravel_index(Teqmask.argmin(), Teqmask.shape)

        dfn = df * 1.0E-6
        doffn = doff * 1.0E6
        # print("Minimum temperature %.3f" %np.amin(Teqmask),"mK occurs at detuning %2.f" %dfn[
        # freqmin],"MHz and offset %.f" %doffn[offsetmin],"micron")

        CS = plt.figure(figsize=(12, 10))
        cmap = plt.cm.Blues_r
        CS = plt.pcolor(dfn, doffn, np.log10(Teqmask), cmap=cmap,
                        vmin=np.log10(mintemp),
                        vmax=np.log10(maxtemp))
        CS = plt.ylim(doffn[0], doffn[-1])
        CS = plt.xlim(dfn[0], dfn[-1])
        CS = plt.xlabel('Laser Detuning $\Delta$ (MHz)')
        CS = plt.ylabel('Laser Offset $d$ $(\mu$ m$)$ ')
        CS = plt.title(title, y=1.03)
        cbar = plt.colorbar()
        cmin = np.log10(mintemp)
        cmax = np.log10(maxtemp)
        ticks = np.linspace(cmin, cmax, 8)
        # print(ticks[-1], cmax, 10 ** cmax, np.max(Teqmask))
        cbar.set_ticks(ticks)
        ticks = 10 ** ticks
        ticks = ["%.2f" % label for label in ticks]

        cbar.ax.set_yticklabels(ticks)

        CS = plt.contour(dfn, doffn, Teqmask, colors="white", levels=contourlevels)
        CS = plt.clabel(CS, inline=1, fontsize=20, fmt='%1.2f')
        if showminimum is True:
            CS = plt.plot(dfn[freqmin], doffn[offsetmin], "*", color="white", markersize=20.0,
                          label="Min. Temp. = %.2f at (%2.f,%2.f)" % (
                              np.amin(Teqmask), dfn[freqmin], doffn[offsetmin]))

            leg = plt.legend(loc=2, fontsize=13, numpoints=1, fancybox=False)
            leg.get_frame().set_alpha(1)
        plt.show()
        return CS

    def plot_torque_result(self, torque, df, doff, temp=None, title="Equilibrium Torque",
                           contourlevels=None, colormidline=.5):

        def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
            '''
            Function to offset the "center" of a colormap. Useful for
            data with a negative min and positive max and you want the
            middle of the colormap's dynamic range to be at zero

            Input
            -----
              cmap : The matplotlib colormap to be altered
              start : Offset from lowest point in the colormap's range.
                  Defaults to 0.0 (no lower ofset). Should be between
                  0.0 and `midpoint`.
              midpoint : The new center of the colormap. Defaults to
                  0.5 (no shift). Should be between 0.0 and 1.0. In
                  general, this should be  1 - vmax/(vmax + abs(vmin))
                  For example if your data range from -15.0 to +5.0 and
                  you want the center of the colormap at 0.0, `midpoint`
                  should be set to  1 - 5/(5 + 15)) or 0.75
              stop : Offset from highets point in the colormap's range.
                  Defaults to 1.0 (no upper ofset). Should be between
                  `midpoint` and 1.0.
            '''
            cdict = {
                'red': [],
                'green': [],
                'blue': [],
                'alpha': []
            }

            # regular index to compute the colors
            reg_index = np.linspace(start, stop, 257)

            # shifted index to match the data
            shift_index = np.hstack([
                np.linspace(0.0, midpoint, 128, endpoint=False),
                np.linspace(midpoint, 1.0, 129, endpoint=True)
            ])

            for ri, si in zip(reg_index, shift_index):
                r, g, b, a = cmap(ri)

                cdict['red'].append((si, r, r))
                cdict['green'].append((si, g, g))
                cdict['blue'].append((si, b, b))
                cdict['alpha'].append((si, a, a))

            newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
            plt.register_cmap(cmap=newcmap)

            return newcmap

        plt.rcParams['font.size'] = 20
        Trq = torque
        Trqmask = np.ma.array(Trq, mask=np.isnan(Trq))

        if temp is not None:
            Teqmask = np.ma.array(temp, mask=np.isnan(temp))

        dfn = df * 1.0E-6
        doffn = doff * 1.0E6

        CS = plt.figure(figsize=(15, 10))
        cmap = plt.cm.RdYlBu
        cmap = shiftedColorMap(cmap, midpoint=colormidline)
        CS = plt.pcolor(dfn, doffn, Trqmask, cmap=cmap, vmin=np.nanmin(Trqmask),
                        vmax=np.nanmax(Trqmask))
        CS = plt.ylim(doffn[0], doffn[-1])
        CS = plt.xlim(dfn[0], dfn[-1])
        CS = plt.xlabel('Laser Detuning $\Delta$ (MHz)')
        CS = plt.ylabel('Laser Offset $d$ $(\mu$ m$)$ ')
        CS = plt.title(title, y=1.03)
        CS = plt.colorbar()
        if temp is not None:
            CS = plt.contour(dfn, doffn, Teqmask, colors='gray', levels=contourlevels)
            CS = plt.clabel(CS, inline=1, fontsize=14, fmt='%1.2f')

        theline = plt.contour(dfn, doffn, Trqmask, colors="black", levels=[0], linewidths=3,
                              fmt="O Torque")

        fmt = {}
        strs = ['Zero Torque']
        for l, s in zip(theline.levels, strs):
            fmt[l] = s
        plt.clabel(theline, theline.levels[::2], inline=True, fmt=fmt, fontsize=18)

        plt.show()

        return

    def generate_y_plot(self, ueq=None, temp=None, detun=None, offset=None,
                        ymin=None, ymax=None,num=300, normalize=True, plot_intensity=True,
                        plot_scattering=True, plot_torque=True, plot_offset=True,
                        show_params=True):

        if ueq is None and temp is None:
            print("No temperature recieved!")
            return False

        if temp is not None:
            ueq = np.sqrt(temp * 2 * self.kB / self.m)

        rp = self.rp
        wy = self.wy
        vk = self.vk
        wr = self.wr

        constants = (self.hbar * self.k) * self.gamma0 * self.s0 * self.sig0 / np.sqrt(np.pi)

        if detun is None:
            detun = self.det
        if offset is None:
            offset = self.d
        if ymax is None:
            ymax = offset + self.wy * 1.85
        if ymin is None:
            ymin = offset - self.wy * 1.85
        delta = 2. * detun * 2 * np.pi / self.gamma0

        yrange = np.linspace(ymin, ymax, num)
        ystep = yrange[1] - yrange[0]

        ploty = yrange * 1E6

        tors = np.ones(num)
        scts = np.ones(num)

        for i in range(num):
            y = yrange[i]
            xmin = -np.sqrt(rp - y ** 2)
            xmax = np.sqrt(rp - y ** 2)

            scts[i] = integ.dblquad(lambda v, x:
                                    self.density(x, y)
                                    * np.exp(-2 * (y - offset) ** 2 / wy ** 2) *
                                    np.exp(-v ** 2) /
                                    (1 + 2 * self.s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2) +
                                     (delta - (wr * y / vk) - (ueq * v / vk)) ** 2),
                                    xmin, xmax,
                                    lambda v: -1 * np.inf, lambda v: np.inf,
                                    epsabs=1.49e-7, epsrel=1.49e-7)[0]
            tors[i] = scts[i] * constants * y * ystep
            scts[i] = scts[i] * constants * ystep / (self.hbar * self.k)

        figr = plt.figure(figsize=(12, 10))

        if not normalize:
            if plot_torque:
                iMax = np.max(np.max(tors), np.abs(np.min(tors)))

            if plot_scattering:
                iMax = np.max(scts)

            if plot_torque and plot_scattering:
                iMax = np.max((np.max(tors), np.abs(np.min(tors))), np.max(scts))
        else:
            iMax = 1

        if plot_intensity is True:
            modwaist = wy * 1.0E6
            intensitylabel = "Intensity of Waist " + "%G $\mu m$" % modwaist
            intens = np.exp(-2 * (yrange - offset) ** 2 / (wy ** 2))
            figr = plt.plot(iMax * intens, ploty,
                            label=intensitylabel, color='black')
        Torqnorm = ""
        Scatnorm = ""
        if normalize is True:
            torsmax = np.max((np.max(tors), np.abs(np.min(tors))))
            Torqnorm = " (Norm. Factor %G)" % torsmax
            tors = tors / torsmax

            sctsmax = np.max(scts)
            Scatnorm = " (Norm. Factor %G)" % sctsmax
            scts = scts / sctsmax

        if plot_torque is True:
            figr = plt.plot(tors, ploty, label="Torque" + Torqnorm, color='red')

        if plot_scattering is True:
            figr = plt.plot(scts, ploty, label="Scattering" + Scatnorm, color="Blue")
        offmax = 1
        if plot_offset is True:
            if plot_torque is True:
                offmax = np.max([np.max(tors), np.abs(np.min(tors))])
            if plot_scattering is True:
                offmax = np.max(scts)
            if plot_scattering is True and plot_torque is True:
                offmax = np.max((np.max(tors), np.max(scts)))
            offlabel = "Offset of " + str(offset * 1E6) + " $\mu m$"
            figr = plt.plot([-offmax, offmax], [1.0E6 * offset, offset * 1.0E6], label=offlabel,
                            color="gray",
                            linestyle="dashed")
        title = list("Spatial Distribution of")
        if plot_intensity: title += list(" Intensity ")

        if str(title[-1]) is " " and plot_scattering:
            title[-1] = ","
        if plot_scattering: title += list(" Scattering ")

        if str(title[-1]) is " " and plot_torque:
            title[-1] = ","
        if plot_torque: title += list(" Torque")
        title = "".join(title)

        figr = plt.title(title)
        figr = plt.ylabel("y ($\mu m$)")
        figr = plt.grid(True)
        figr = plt.ylim((ymin * 1E6, ymax * 1E6))
        figr = plt.legend(loc=3, fontsize=12)
        figr = plt.plot([-iMax, iMax], [0, 0],
                        color="gray")
        if show_params is True:
            figr = plt.text(-.7 * iMax, 1E6 * ymax * .80, "$\Delta$=%G MHz" % (1.0E-6 * detun),
                            fontsize=20)
            figr = plt.text(-.7 * iMax, 1E6 * ymax * .60, "$d=$ %G $\mu$m" % (offset * 1.0E6),
                            fontsize=20)
            figr = plt.text(-.7 * iMax, 1E6 * ymax * .40,
                            "$T=$ %G mK" % (1000 * ueq ** 2 * self.m / (2 * self.kB)),
                            fontsize=20)

        plt.show()
        return 1


if __name__ == '__main__':
    a = ItanoAnalysis(quiet=False, ywidth=30E-6, saturation=.5)  # ,saturation=.1)

    df = np.linspace(-150.0E6, -1E6, 3)
    doff = np.linspace(0, 20.0E-6, 3)
    A = a.scan_detuning_and_offset(detuninglist=df, offsetlist=doff, no_rotating_wall=True)

