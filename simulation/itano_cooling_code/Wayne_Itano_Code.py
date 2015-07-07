__author__ = 'sbt'

""" Began life as MATLAB code written by Charles Xu.
    Written for Python by Steven B. Torrisi, Summer 2015.
    Contains a class called ItanoAnalysis, based primarily on the PRA Paper by W. Itano from 1988
    "Perpendicular laser cooling of a rotating ion plasma in a Penning trap"
    with modifications which take into account a rotating wall potential, as well as
    some recoil heating from a parallel laser beam frmo the PRA paper by W. Itano from 1982
    "Laser Cooling of Atoms Stored in Harmonic and Penning Traps".

    July 6 Revision

"""
import numpy as np
import scipy.optimize as opt
import scipy.integrate as integ
import matplotlib.pyplot as plt
from scipy.constants import pi

plt.rcParams['font.size'] = 16


class ItanoAnalysis:
    """
    Contains the tools to conduct calculations of total energy balance, torque, scattering events
    from a perpendicular doppler cooling beam acting incident on a Berylium ion plasma.

    Methods inside the function allow one to find, for a plasma of a given radius and density,
    the thermal equilibrium which results from a perpendicular doppler cooling laser acting on the
    plasma.
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
    pardiff = 2 * np.pi * 9E6

    def __init__(self,
                 defaultoff=30.0E-6, defaultdet=-500E6, wr=2 * pi * 45.0E3, Tguess=1E-3,
                 saturation=.5,
                 dens=2.77E9, ywidth=2.0E-6, radius=225.0E-6, quiet=True, spar=.2, wpar=9E6):

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
        self.u0 = np.sqrt(
            2 * self.kB * Tguess / self.m)  # initial thermal velocity used in Boltzmann factor
        self.s0 = saturation  # saturation parameter for laser interaction (set to 0 to turn off saturation effects)
        self.sig0 = dens
        self.wy = ywidth


        # if you choose to include a parallel laser

        self.spar = spar
        self.wpar = wpar

        self.counter = 0

    def density(self, x, y, radius=None):
        """
        Defines a density scaling which returns 0 outside of the circular bounds of the plasma.
        Useful because it allows for rectangular integration bounds (i.e. [-r,r])
        instead of complicated functional bounds which are difficult to work with
        (such as [-\sqrt{1+x^2}, \sqrt{1+x*2}] for example).

        uses the plasma radius


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
                            ((v) + 5 * self.rnorm / (3 * u)) * np.exp(-v ** 2) /
                            ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                              + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                            -np.inf, np.inf,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * rp, lambda x, y: rp)
        self.counter += 1
        return ret[0]

    def dEavgAndParallel(self, u, detun=None, offset=None):
        """

        Returns the average rate of change of energy from the perpendicular doppler cooling laser,
        including the effects of a laser in the z-direction (thus including a parallel cooling laser
        too. This contribution comes in the form of a constant recoil heating in the x and y degrees of freedom).

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

        alpha = (2 / 9) * self.hbar * self.k * rp ** 2 * np.pi ** (3 / 2) * self.spar \
                / (s0 * self.m * (1 + 2 * self.spar + (2 / self.gamma0) ** 2) * self.pardiff ** 2)

        ret = integ.tplquad(lambda y, x, v:
                            self.density(x, y)
                            * np.exp(-(y - offset) ** 2 / wy ** 2) *
                            ((v) + 5 * self.rnorm / (3 * u)) * np.exp(-v ** 2) /
                            ((1 + 2 * s0 * np.exp(-2 * (y - offset) ** 2 / wy ** 2)
                              + (delta - (wr * y / vk) - (u * v / vk)) ** 2)),
                            -np.inf, np.inf,
                            lambda x: -1 * rp, lambda x: rp,
                            lambda x, y: -1 * rp, lambda x, y: rp)

        # print("Integral Evaluated",ret[0],"with temperature",u**2*8.96*1.673E-27/(2*1.38E-23))
        self.counter += 1
        return u * ret[0] + alpha

    def totalscatterperp(self, ueq, detun=None, offset=None):
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
                                 get_temp=True, get_torque=False, get_total_scatter=False,
                                 plot=False,
                                 detuninglist=None, offsetlist=None, parheating=False) -> list:
        """
            Conducts a two-dimensional scan across the detuning/offset
            parameter space.
            Returns an array of output elements.

            if plot is true, though this feature is not in yet, it can
            produce a nice plot for you.


            :rtype : list which contains 3 elements:
                1) Teq element, a 2d array of the equilibrium temperatures,
                    or the integer 0 if get_temp is false.
                2) Trq element, a 2d array of the net torques, or the integer
                    0 if get_torque is false.
                3) Sct element, a 2d array of the net scattering rate, or the
                    integer 0 if get_total_scatter is false.

            :param detmin: Minimum detuning. Make negative, units of Hz.
            :param detmax: Maximum detuning. Make negative, units of Hz.
            :param detN:  Number of detunings to evaluate
                            (interpolates linearly)
            :param offmin: Minimum offset. Make nonnegative, units of meters.
            :param offmax: Maximum offset. Make nonnegative, units of meters.
            :param offN:  Number of offsets to evaluate (interpolates linearly)
            :param get_temp: If true, returns temperature in the result array.
            :param get_torque: If true, returns the torques in the result array.
            :param get_total_scatter: If true, returns the scattering rates.
            :param detuninglist: Allows for the use of a custom list of
                            detunings instead of a linear interpolation.
            :param offsetlist: Allows for the use of a custom list of offsets instead
                                        of a linear interpolation.
            :param plot:
            """
        Tmax = 1.0E-1
        Tmin = 1.0E-5
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        Trq = 0
        Sct = 0

        if parheating is False:
            dEfunction =self.dEavg
        else:
            dEfunction = self.dEavgAndParallel

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
        if get_torque is True:
            Trq = np.zeros([len(doff), len(df)])
        if get_total_scatter is True:
            Sct = np.zeros([len(doff), len(df)])
        for W in range(len(df)):
            for D in range(len(doff)):
                if self.quiet is False:
                    print("-------------")
                    print("Feeding in detuning/offset", dw[W] / (2 * 3.14159),
                          doff[D], W * (len(doff)) + D + 1, "of",
                          len(df) * len(doff), "Evaluations",
                          100 * (W * (len(doff)) + D + 1) / (len(df) * len(doff)), "%")

                # The difference between each of these three tries is that
                try:
                    U[D][W] = opt.brentq(dEfunction, umin, umax,
                                         args=(dw[W], doff[D]), xtol=1e-4, rtol=3.0e-7)
                except:
                    try:
                        U[D][W] = opt.brentq(dEfunction, umin * .1, umax * 10,
                                             args=(dw[W], doff[D]), xtol=1e-4,
                                             rtol=3.0e-7)
                    except:

                        try:
                            U[D][W] = opt.brentq(dEfunction, umin * .01, umax * 100,
                                                 args=(dw[W], doff[D]), xtol=1e-4,
                                                 rtol=3.0e-7)
                        except:
                            try:
                                U[D][W] = opt.brentq(dEfunction, umin * .001, umax * 1000,
                                                     args=(dw[W], doff[D]), xtol=1e-5,
                                                     rtol=3.0e-7)

                            except:
                                U[D][W] = np.nan  # Obviously an error occured
                                print(
                                    "Warning! The integral was not able to evaluate for some reason. "
                                    "Storing u as nan and moving on...")

                if get_torque:
                    Trq[D][W] = self.totaltorque(U[D][W], dw[W], doff[D])

                if get_total_scatter:
                    Sct[D][W] = self.totalscatterperp(U[D][W], dw[W], doff[D])

                if self.quiet is False:
                    print("u/T recorded:", U[D][W],
                          ((U[D][W]) * U[D][W] * 8.96 * 1.673E-27 / (2 * 1.38E-23)), "after",
                          self.counter, "integral calls.")

                self.counter = 0
        if get_temp is True:
            Teq = np.array([u ** 2 * self.m / (2 * self.kB) for u in U])
        else:
            Teq = 0
        return [Teq, Trq, Sct]

    def getueq(self, detun=None, offset=None, returntemperature=False):
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

        Tmax = 1.0E-1
        Tmin = 1.0E-5
        umax = np.sqrt(Tmax * self.kB * 2 / self.m)
        umin = np.sqrt(Tmin * self.kB * 2 / self.m)
        ret = opt.brentq(self.dEavg, umin, umax, args=(detun, offset), xtol=1e-4, rtol=3.0e-7)
        if returntemperature is False:
            return ret[0]
        else:
            return ret[0] ** 2 * self.m / (2 * self.kB)

    def plot_temperature_result(self, temp, df, doff, showminimum=True,
                                title="Equilibrium Temperature",contourlevels=None,
                                maxtemp=None):
        """
        Input the 0th array which returns from scan detuning and offset to
        make a plot of the temperature.

        Assumptions about input parameters are as follows:


        """

        plt.rcParams['font.size'] = 16
        Teq = temp[0] * 1000

        Teqmask = np.ma.array(Teq, mask=np.isnan(Teq))
        if maxtemp is not None:
            Teqmask = np.ma.array(Teqmask, mask=Teq >= maxtemp)

        offsetmin, freqmin = np.unravel_index(Teqmask.argmin(), Teqmask.shape)

        dfn = df * 1.0E-6
        doffn = doff * 1.0E6

        # print("Minimum temperature %.3f" % np.amin(Teqmask),
        #      "mK occurs at detuning %2.f" % dfn[freqmin], "MHz and offset %.f" % doffn[offsetmin],
        #      "micron")

        CS = plt.figure(figsize=(12, 10))
        cmap = plt.cm.Blues_r
        CS = plt.pcolor(dfn, doffn, Teqmask, cmap=cmap, vmin=np.nanmin(Teqmask),
                        vmax=np.nanmax(Teqmask))
        CS = plt.ylim(doffn[0], doffn[-1])
        CS = plt.xlim(dfn[0], dfn[-1])
        CS = plt.xlabel('Laser Detuning (MHz)')
        CS = plt.ylabel('Laser Offset $(\mu m)$ ')
        CS = plt.title(title,
                       y=1.03)
        CS = plt.colorbar()

        if contourlevels is not None:
            CS = plt.contour(dfn, doffn, Teqmask, colors='white', levels=contourlevels)
            CS = plt.clabel(CS, inline=1, fontsize=10, fmt='%1.2f')

        if showminimum is True:
            CS = plt.plot(dfn[freqmin], doffn[offsetmin], "*", color="white", markersize=20.0,
                          label="Min. Temp. = %.3f at (%2.f,%2.f)" % (
                              np.amin(Teqmask), dfn[freqmin], doffn[offsetmin]))
            leg = plt.legend(loc=2, fontsize=13, numpoints=1, fancybox=False)
            leg.get_frame().set_alpha(1)
        CS=plt.show()
        return CS



if __name__ == '__main__':
    a = ItanoAnalysis(quiet=False, spar=1.0E30)

    A = a.scan_detuning_and_offset(detN=2, offN=2, get_torque=True, get_total_scatter=True,
                                   parheating=False)
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
