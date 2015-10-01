import sys
import time
import subprocess
from contracts import contract

from scipy.constants import pi
import numpy as np
import scipy.constants as u
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import numpy.matlib as matlib
import matplotlib.gridspec as gridspec


"""
Contains the ModeAnalysis class, which can simulate the positions of ions in a crystal
of desired size. The class contains methods for the generation of a crystal,
relaxation to a minimum potential energy state, and determination of axial and (eventually) planar modes of motion
by methods derived by Wang, Keith, and Freericks in 2013.

Translated from MATLAB code written by Adam Keith by Justin Bohnet.
Standardized and slightly revised by Steven Torrisi.

Be careful. Sometimes when using the exact same parameters this 
code will make different crystals with the same potential energy. That is,
crystal are degenerate when reflecting over the axes.
"""


class HexLattice:
    """
    Class defining methods relating to a perfect, closed 2D hexagonal lattice.
    """
    def __init__(self, s: 'int,>=0', scale=1):
        """

        :param s: num of closed shells (int)
        :param scale: target vertex-vertex spacing (float)
        :return: none
        """
        self._s = s
        self._scale = scale

    @property
    @contract(returns='tuple(array[N],array[N]),N>0')
    def xy(self):
        """coordinates of lattice vertices

        :return: (x, y)
        """
        self._x, self._y = HexLattice.hex_lattice(self._s, self._scale)
        return self._x, self._y

    @staticmethod
    @contract(shells='int,>=0', scale='>0', returns='tuple(array[N],array[N]),N>0')
    def hex_lattice(shells, scale=20e-5):
        """Generate closed shell hexagonal lattice with shells and scale spacing.

        :param scale: scales lattice
        :return: x, y coordinates of a hexagonal lattice
        """
        nion = 1 + 6 * np.sum(range(1, shells + 1))
        shellsx = []
        shellsy = []
        for s in range(0, shells + 1):
            x, y = HexLattice.hex_shell(s)
            shellsx.append(x)
            shellsy.append(y)
        return np.hstack(shellsx)*scale, np.hstack(shellsy)*scale

    @staticmethod
    @contract(s='int,>=0', returns='tuple(array[N],array[N]),N>0')
    def hex_shell(s: int):
        """ Generates the sth shell of a 2-d hexagonal lattice.

        :param s: shell number

        :return: (x, y) vertex coordinates
        """
        if s==0:
            x = np.array([0])
            y = np.array([0])
        else:
            a = list(range(s, -s - 1, -1))
            a.extend(-s * np.ones(s - 1))
            a.extend(list(range(-s, s + 1)))
            a.extend(s * np.ones(s - 1))

            b = list(range(0, s + 1))
            b.extend(s * np.ones(s - 1))
            b.extend(list(range(s, -s - 1, -1)))
            b.extend(-s * np.ones(s - 1))
            b.extend(list(range(-s, 0)))

            x = np.sqrt(3) / 2.0 * np.array(b)
            y = 0.5 * np.array(b) + np.array(a)
        return x, y

    @staticmethod
    @contract(nshells='int,>=0', returns='>=0')
    def get_nvert_from_nshells(nshells):
        """Number of vertices for a fully filled hexagonal lattice
        with nshells. A single lattice site lies at the center.

        :param nshells: number of hexagonal shells
        :return: number of vertices
        """
        return 1 + 6 * np.sum(range(1, nshells + 1))

    @staticmethod
    @contract(nvert='int,>=0', returns='>=0')
    def get_nshells_from_nvert(nvert):
        """Minimum number of closed shells for hexagonal lattice
        that contains at least nvert vertices.

        :param nvert: number of vertices
        :return: number of shells
        """
        for s in range(0, 20):
            ntest = HexLattice.get_nvert_from_nshells(s)
            if nvert < ntest:
                return s-1

    def plot(self):
        plt.clf()
        for si in range(self._s):
            x, y = self.hex_shell(si)
            x = np.hstack((x, x[0]))
            y = np.hstack((y, y[0]))
            plt.plot(x, y, 'k.', markersize=30)
            plt.plot(x, y, 'b-')
            plt.axis('equal')
        plt.axis('off')
        plt.show()

class PenningTrap:
    """
    Penning trap configuration and potentials.
    """
    _electrodes_applied_v = [0, -1750, -2000]
    # see Carson Teale's final paper
                   # [endcap,    mid,       center]
    _eerz_2013teal = np.array(\
                    [[0.0756,   0.5157,    0.4087],      # C0
                    [-0.0001,  -0.005,     0.005],       # C1
                    [ 1.9197e3, 3.7467e3, -5.6663e3],    # C2
                    [ 0       , 0       ,  0       ],    # C3
                    [ 0.6738e7,-5.3148e7,  4.641e7]])    # C4
    _eerphi_wall_2013teale = {2: 1612, 3: 30010} # index is wall_order
    def __init__(self, b:'>0'=4.4588,
            electrodes_applied_v: 'array[1x3]'=_electrodes_applied_v,
            electrodes_expansion_in_rz: 'array[3x4]'=_eerz_2013teal,
            wall_f: '>0'=180e3, wall_v: '>=0'=5,
            wall_r: '>=0'=0.01, wall_order: 'int,>0,<4'=2):
        """
        :param b: magnetic field along z-axis (Tesla)
        :param electrodes_applied_v: [end cap, middle, center]
                potential on electrodes (Volts)
        :param electrodes_expansion_in_rz: array containing
                expansion (in cylindrical
                coordinates) about z=0,r=0 of potential due to 1 V on each
                electrode
        :param wall_f: rotating wall frequency (Hz)
        :param wall_v: rotating wall potential (V)
        :param wall_r: rotating wall distance to trap center (m)
        :param wall_order: rotating wall order {1,2,3}
        :param verbose: (bool)
        :return: none
        """
        self._b = b
        self._electrodes_expansion_in_rz = \
            np.array(electrodes_expansion_in_rz)
        self._electrodes_applied_v = np.array(electrodes_applied_v)

        self._wall_v = wall_v
        self._wall_f = wall_f
        self._wall_order = wall_order
        self._wall_r = wall_r

    def __str__(self):
        s = "b={} T\n".format(self.b)
        s+= "electrodes_electrodes_expansion_in_rz={} V [end, mid, center]\n"\
            .format(self.trap_v)
        s+= "electrodes_applied_v={}\n".format(self.electrodes_applied_v)
        s+= "wall_f={} kHz  wall_v={} V  wall_order={}  wall_r={} m"\
            .format(self.wall_f/1e3, self.wall_v, self.wall_order, self.wall_r)
        return s

    @property
    def b(self) -> '>0':
        return self._b

    @property
    def wall_v(self) -> '>0':
        return self._wall_v

    @property
    def wall_f(self) -> '>0':
        return self._wall_f

    @property
    def wall_order(self) -> 'int,>0,<4':
        return self._wall_order

    @property
    def wall_r(self) -> '>0':
        return self._wall_r

    @property
    def electrodes_expansion_in_rz(self):
        """array containing expansion (in cylindrical
        coordinates) about z=0,r=0 of potential due to 1 V on each
        electrode
        """
        return self._electrodes_expansion_in_rz

    @property
    def electrodes_applied_v(self) -> 'array[3]':
        return self._electrodes_applied_v

    @property
    def trap_pot_expansion_in_rz_coef(self):
        """Expansion (in cylindrical coordinates) about z=0, r=0 of
        potential at trap center. This is a sum over all
        electrodes assuming self.trap_v is applied.
        """
        return np.dot(self.electrodes_expansion_in_rz,
                      self.electrodes_applied_v)  # formerly Coeff[]

    def trap_pot_at_zr(self, z, r):
        """
        Trap electronic potential at (r,z). Does not include rotating wall.
        """
        # rpowpleg is terms of r**n * Legendre(n, z/r)
        rpowpleg = np.array([1,
                              z,
                              0.5*(-r**2+3*z**2),
                              0.5*(-3*r**2*z+5*z**3),
                              0.125*(3*r**4-30*r**2*z**2+35*z**4)])
        v = np.dot(self.trap_pot_expansion_in_rz_coef, rpowpleg)
        return v

    def wall_pot_expansion_in_rphi_coef(self, wall_order,
                eerphi=_eerphi_wall_2013teale):
        """Rotating wall potential can be expanded near the trap center
        in terms of r, phi. Carson Teal obtained an approximate form for this
        in 2013.  See Teale final report Aug 6, 2013.
        Note that this assumes z=0.
        """
        if wall_order==2 or wall_order==3:
            return eerphi[wall_order]
        else:
            # unsupported wall order
            sys.exit(-225)

    def wall_pot_at_zrphi(self, z, r, phi, wall_order, method='teale2013'):
        """
        Rotating wall electronic potential at (z,r,phi).

        :param z: z coordinate (meters)
        :param r: radial coordinate (meters)
        :param phi: azimuthal angle (turns, 0:2*np.pi)
        """
        # The most complete representation of the wall potential would be to
        # a) apply 1 V to one rotating wall electrode
        # b) use BEM to compute the potential in volume near trap center
        # c1) save these potentials for fully-general calculations later
        # c2) in z=0 plane, use expansion V(r, phi) = sum A(l,m) r^l cos(m phi)
        #     see John's note 2/19/2015 in Onenote
        # c3) whatever Carson did
        if method == 'full':
            # this is c1 -- needs to be done!
            sys.exit(-222)
        elif method == 'bem_inplane':
            # this is c2 -- needs to be done!
            if z!= 0:
                sys.exit(-223)
            else:
                sys.exit(-221)
        elif method == 'teale2013':
            alm = self.wall_pot_expansion_in_rphi_coef(wall_order)
            x = r*np.cos(phi)
            y = r*np.sin(phi)
            if wall_order==2:
                pot = alm*(x**2 - y**2)
            if wall_order==3:
                pot = (x**3 - 3*x*y**2)
        else:
            sys.exit(-220)
        return pot

    def plot_wall_pot_rphi(self, method='teale2013'):
        r = np.linspace(0, 2, 50)*1e-3
        phi = np.linspace(0, 2*np.pi, 360)
        mr, mphi = np.meshgrid(r, phi)
        ax1 = plt.subplot(1, 2, 1, projection="polar")
        pot = self.wall_pot_at_zrphi(z=0, r=mr, phi=mphi,
                                     wall_order=2, method=method)
        ax1.pcolormesh(mphi, mr, pot, cmap=plt.cm.bwr)
        ax2 = plt.subplot(1, 2, 2, projection="polar")
        pot = self.wall_pot_at_zrphi(z=0, r=mr, phi=mphi,
                                     wall_order=3, method=method)
        ax2.pcolormesh(mphi, mr, pot, cmap=plt.cm.bwr)
        plt.show()

    def plot_trap_pot_along_z(self):
        # extract trap potential, pot
        zs = np.linspace(-20, 20, 50)*1e-3
        trap_pot = [g.trap.trap_pot_at_zr(z, 0) for z in zs]
        # calculate potential for perfectly quadratic trap
        c0 = self.trap_pot_expansion_in_rz_coef[0]
        c2 = self.trap_pot_expansion_in_rz_coef[2]
        r = 0
        quad_pot = c0 + [c2*0.5*(-r**2+3*z**2) for z in zs]
        # plot
        plt.plot(zs, trap_pot)
        plt.plot(zs, quad_pot)

    def plot_trap_pot_along_r(self):
        # extract trap potential, pot
        rs = np.linspace(-20, 20, 50)*1e-3
        z = 0
        trap_pot = [g.trap.trap_pot_at_zr(z, r) for r in rs]
        # calculate potential for perfectly quadratic trap
        c0 = self.trap_pot_expansion_in_rz_coef[0]
        c2 = self.trap_pot_expansion_in_rz_coef[2]

        quad_pot = c0 + [c2*0.5*(-r**2+3*z**2) for r in rs]
        # plot
        plt.plot(rs, trap_pot)
        plt.plot(rs, quad_pot)


class Geonium():
    """
    Single ion confined in a Penning trap.
    """
    _amu = 1.66057e-27
    _q = 1.602176565E-19
    _m_Be = 9.012182 * _amu
    _k_e = 8.9875517873681764E9

    def __init__(self, trap: PenningTrap,
                          m: '>0'=_m_Be, q: '>0'=_q):
        """Physics related to a single isolated ion in Penning trap

        :param trap: instance of object PenningTrap
        :param m: ion mass (kilogram)
        :param q: ion charge (coulomb)
        :return: none
        """
        self.trap = trap
        self._m = m
        self._q = q
        c = self.trap.trap_pot_expansion_in_rz_coef
        self._ax_f = np.sqrt(2*q*c[2]/m)/2/np.pi
        self._cyc_f = q*self.trap.b/(2*np.pi*m)
        self._mag_f = \
            0.5*(self.cyc_f - np.sqrt(self.cyc_f**2 - 2*self.ax_f ** 2))/2/np.pi

        self.check_magnetron()

        # quadratic voltage at trap center
        self.V0 = (0.5*m*self.ax_f**2)/q/2/np.pi

    @property
    def m(self):
        return self._m

    @property
    def q(self):
        return self._q

    @property
    def cyc_f(self):
        return self._cyc_f

    @property
    def ax_f(self):
        return self._ax_f

    @property
    def mag_f(self):
        return self._mag_f

    def __str__(self):
        s = "m={:.2e} kg  q={:.2e}\n".format(self.m, self.q)
        s+= "cyc_f={:.1f} MHz  ax_f={:.2f} MHz  mag_f={:.1f} kHz"\
            .format(self.cyc_f/1e6, self.ax_f/1e6, self.mag_f/1e3)
        return s

    def check_magnetron(self):
        """Penning does not trap for ion rotation frequency below magnetron freq.

        :param ions: ModeAnalysis object
        :return: bool
        """
        if self.mag_f > self.trap.wall_f:
            print("fmag ({:.0f}) > fwall ({:.0f})"
                  .format(self.mag_f, self.trap.wall_f))
            return False
        else:
            return True

    def __theory_unit_converter(self):
        _length = ((self._k_e * self.q**2)/(.5 * self.m * (self.ax_f*2*np.pi)**2))**(1/3)
        _time = 1/(2*np.pi*self.ax_f)
        _velocity = _length/_time
        _energy = .5*self.m*(self.ax_f*2*np.pi)**2*_length**2
        _wall = self.trap.wall_f/self.ax_f
        _cyclotron = self.cyc_f/self.ax_f
        _mass = 1
        lookup = {
            "length": _length,      # characteristic length
            "time": _time,          # characteristic time
            "velocity": _velocity,  # characteristic velocity
            "energy": _energy,      # characteristic energy
            "wall": _wall,          # dimensionless rotation
            "cyclotron": _cyclotron,  # dimensionless cyclotron
            "mass": _mass           # dimensionless mass
            }
        return lookup

    def to_theory_units(self, unit, value):
        lookup = self.__theory_unit_converter()
        return value/lookup[unit]

    def to_si_units(self, unit, value):
        lookup = self.__theory_unit_converter()
        return value*lookup[unit]

    def __pot_energy(self, x: 'array[M],M>0', y: 'array[M],M>0') -> 'float':
        """
        Computes the potential energy of the single ion, taking into consideration:
        Coulomb repulsion
        qv x B forces
        Trapping potential
        and some other things (#todo to be fully analyzed; june 10 2015)

        :param x: ion coordinate (meters)
        :param y: ion coordinate (meters)
        :return: potential energy (volts)
        """
        # Frequency of rotation, mass and the number of ions in the array
        wr = 2*np.pi*self.wall_f
        m = self.m
        q = self.q
        b = self.b
        k_e = self.c_k_e
        wa2 = self.wall_coef2
        wa3 = self.wall_coef3
        mult = 1e14  # mystery multiplicative factor

        # One half times the rotational force, the charge times the coeff,
        potential = 0.5*(-m*wr**2 - q*self.pot[2] + q*b*wr) * np.sum((x**2 + y**2)) \
            - q*self.pot[3] * np.sum((x**2 + y**2)**2) \
            + np.sum(wa2*(x**2 - y**2)) \
            + np.sum(wa3*(x**3 - 3*x*y**2))
        return mult*potential


class IonCrystal2d():
    _crystal_scale_fudge_factor = 20e-6

    def __init__(self, geonium: Geonium,
                     nion: 'int,>0'=27, perturb=True):
        """2-d array of ions confined in a Penning trap
        Assume all ions are the same mass

        :param geonium: instance of class Geonium
        :param nion: number of ions
        :param perturb: stochastic perturbation during crystal solution (bool)
        :return: none
        """
        self.geonium = geonium
        self._nion = nion
        self._perturb = perturb

        self._x = np.empty(nion)
        self._y = np.empty(nion)
        self._x, self._y = self.__solve_for_minimum_energy_crystal()

    @property
    def nion(self):
        return self._nion

    @property
    def xy(self) -> 'tuple(array[N],array[N]),N>0':
        """Ion coordinates in lab frame

        :return: (x, y) (meters)
        """
        return self._x, self._y

    @property
    def dij(self) -> 'tuple(array[N],array[N]),N>0':
        """ion-ion separation

        :return: dij (NxN)
        """
        x, y = self.xy
        nion = self.nion
        xrep = matlib.repmat(x, nion, 1)
        yrep = matlib.repmat(y, nion, 1)
        dij = np.sqrt((xrep-xrep.T)**2 + (yrep-yrep.T)**2)
        return dij

    @property
    def dij_min(self) -> 'float':
        dij = self.dij
        nion = self.nion
        # transform Jij and dij to 1D vectors for plotting
        dij_row = np.squeeze(np.array(np.reshape(dij, (1, nion**2))))
        dij_nonzero = [d if d>1e-9 else np.NaN for d in dij_row]
        return np.nanmin(dij_nonzero)

    @property
    def dij_max(self) -> 'float':
        dij = self.dij
        nion = self.nion
        # transform Jij and dij to 1D vectors for plotting
        dij_row = np.squeeze(np.array(np.reshape(dij, (1, nion**2))))
        dij_nonzero = [d if d>1e-9 else np.NaN for d in dij_row]
        return np.nanmax(dij_nonzero)

    def __str__(self):
        s = "nion={}  dij_min={:.1f} um  dij_max={:.1f} um"\
            .format(self.nion, self.dij_min*1e6, self.dij_max*1e6)
        return s

    def __solve_for_minimum_energy_crystal(self, perturb=False, timeout=1)\
            -> 'tuple(array[N],array[N]),N>0':
        """
        :param perturb: re-solve for equilibrium after perturbing ions
        :param timeout: time permitted for solving (seconds)
        :return: none
        """
        t00 = time.time()
        # is ion number compatible with filled lattice?
        n_hex_shells = HexLattice.get_nshells_from_nvert(self.nion)
        if HexLattice.get_nvert_from_nshells(n_hex_shells) == self.nion:
            lattice00 = HexLattice(n_hex_shells,
                                  scale=self._crystal_scale_fudge_factor)
            x00, y00 = lattice00.xy
            # initial potential energy
            pot00 = self.__pot_energy(np.hstack((x00,y00)))
        else:
            x00, y00 = self.__generate_lattice_for_arbitrary_n(self.nion)

        # first round of minimization
        x, y = self.__pot_minimize(x00, y00)

        # try multiple rounds of minimization, giving random nudges
        if perturb:
            pot0 = self.__pot_energy(x, y)
            x0 = x
            y0 = y
            strength = 0.01
            epsilon = 1e-6
            dpot = 1e10
            while time.time() - t00 < timeout:
                px = np.dot(strength*np.random.randn(self.nion), x)
                py = np.dot(strength*np.random.randn(self.nion), y)
                x, y = self.__pot_minimize(px)
                pot = self.__pot_energy(x, y)
                if pot < pot0:
                    x0 = x
                    y0 = y
                    dpot = (pot0-pot)/pot00
                    pot0 = pot
                if dpot < epsilon:
                    break
        return x, y

    def __generate_lattice_for_arbitrary_n(self, nion: '>0')\
            -> 'tuple(array[N],array[N]),N>0':
        """Generate lattice for an arbitrary number of ions.
        Homogeneous mass.

        :param nion: number of ions in lattice
        :return: x, y coordinates of vertices
        """
        # number of closed shells
        s = int((np.sqrt(9 - 12 * (1 - nion)) - 3) / 6)
        lattice0 = HexLattice(s)
        x0, y0 = lattice0.xy
        n0 = x0.size
        n_add = nion - n0  # Number of ions left to add
        self._nion = n0

        pair = HexLattice.hex_shell(s + 1)  # generate next complete shell
        xadd = pair[0::2]
        yadd = pair[1::2]

        for i in range(n_add):
            # reset number of ions to do this calculation
            nion += 1

            # make masses all one (add defects later)
            self.md = np.ones(nion)

            v = []  # list to store potential energies from calculation

            # for each ion left to add, calculate potential energy if that
            # ion is added
            for j in range(len(xadd)):
                v.append(self.__pot_energy(
                    np.hstack((x0, xadd[j])),
                    np.hstack((y0, yadd[j])))
                )
            ind = np.argmin(v)  # ion added with lowest increase in potential

            # permanently add to existing crystal
            x0 = np.append(x0, xadd[ind])
            y0 = np.append(y0, yadd[ind])

            # remove ion from list to add
            xadd = np.delete(xadd, ind)
            yadd = np.delete(yadd, ind)

        return x0, y0

    def __pot_minimize(self, x0: 'array[N],N>0', y0: 'array[N],N>0',
                       method="bfgs") -> 'tuple(array[N],array[N]),N>0':
        """Minimize ion crystal potential energy

        :param x0: The position vector which defines the crystal.
        :param y:
        :param method: {"bfgs", "newton"}
        :return: The equilibrium position vector.
        """
        newton_tolerance = 1e-34
        bfgs_tolerance = 1e-34
        nion = self.nion
        u0 = np.hstack((x0, y0))

        if method is "newton":
            result = optimize.minimize(self.__pot_energy, u0, method='Newton-CG',
                jac=self.__force, hess=self.__hessian,
                options={'xtol': newton_tolerance, 'disp': not self.quiet})
        elif method is "bfgs":
            result = optimize.minimize(self.__pot_energy, u0, method='BFGS',
                jac=self.__force,
                options={'gtol': bfgs_tolerance, 'disp': False})
        else:
            sys.exit(-2)
        pos_array = result.x
        return pos_array[:nion], pos_array[nion:]

    def __pot_energy(self, pos_array):
        """
        Computes the potential energy of a 2-d plane of ions,
        taking into consideration:
            Coulomb repulsion
            qv x B forces
            Trapping potential
            and some other things (#todo to be fully analyzed; june 10 2015)

        :param x: ion coordinates (meters)
        :param y: ion coordinates (meters)
        :return: potential energy (volts)
        """
        nion = self.nion
        x = pos_array[0:nion]
        y = pos_array[nion:]

        dx = x.reshape((x.size, 1)) - x  # row vector
        dy = y.reshape((y.size, 1)) - y  # row vector
        dij = np.sqrt(dx ** 2 + dy ** 2)  # distance between ions
        wr = self.geonium.to_theory_units('wall', self.geonium.trap.wall_f)
        wc = self.geonium.to_theory_units('cyclotron', self.geonium.cyc_f)
        m = self.geonium.to_theory_units('mass', self.geonium.m)
        q = self.geonium.q
        wall_order = self.geonium.trap.wall_order
        wc2 = (wall_order==2)*q*self.geonium.trap.wall_pot_expansion_in_rphi_coef(2)
        wc3 = (wall_order==3)*q*self.geonium.trap.wall_pot_expansion_in_rphi_coef(3)
        b = self.geonium.trap.b
        k_e = self.geonium._k_e
        # TODO: why aren't x and y converted to theory units?

        with np.errstate(divide='ignore'):
            vc = np.where(dij != 0., 1 / dij, 0)

        pot = 0.5*(-m*wr**2 - q*wc2 + q*b*wr) * np.sum((x**2 + y**2)) \
            - q*wc3 * np.sum((x**2 + y**2)**2) \
            + np.sum(wc2 * (x**2 - y**2)) \
            + np.sum(wc3 * (x**3 - 3*x*y**2)) \
            + 0.5 * k_e*q**2 * np.sum(vc)
        return pot

    def __hessian(self, pos_array):
        """Calculate Hessian of potential"""

        x = pos_array[0:self.n]
        y = pos_array[self.n:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx ** 2 + dy ** 2)
        with np.errstate(divide='ignore'):
            rsep5 = np.where(rsep != 0., rsep ** (-5), 0)
        dxsq = dx ** 2
        dysq = dy ** 2

        wr = self.to_theory_units('wall', self.wall_f)
        wc = self.to_theory_units('cyclotron', self.cyc_f)
        m = self.to_theory_units('mass', self.m)
        q = self.q
        wall_order = self.geonium.trap.wall_order
        wc2 = (wall_order==2)*q*self.geonium.trap.wall_pot_expansion_in_rphi_coef(2)
        wc3 = (wall_order==3)*q*self.geonium.trap.wall_pot_expansion_in_rphi_coef(3)
        b = self.b
        k_e = self._k_e

        # X derivatives, Y derivatives for alpha != beta
        hess_xx = np.mat((rsep ** 2 - 3 * dxsq) * rsep5)
        hess_yy = np.mat((rsep ** 2 - 3 * dysq) * rsep5)

        # Above, for alpha == beta
        hess_xx += np.mat(np.diag(-2*m*(wr**2 - wr*wc + .5 - wc2) -
                              np.sum((rsep**2 - 3*dxsq)*rsep5, axis=0)))
        hess_yy += np.mat(np.diag(-2*m*(wr**2 - wr*wc + .5 + wc2) -
                              np.sum((rsep**2 - 3*dxsq)*rsep5, axis=0)))

        # Mixed derivatives
        hess_xy = np.mat(-3*dx*dy*rsep5)
        hess_xy += np.mat(np.diag(3 * np.sum(dx*dy*rsep5, axis=0)))

        hess = np.bmat([[hess_xx, hess_xy], [hess_xy, hess_yy]])
        hess = np.asarray(hess)
        return hess

    def __force(self, pos_array):
        """
        Computes the net forces acting on each ion in the crystal;
        used as the jacobian by find_eq_pos to minimize the potential energy
        of a crystal configuration.

        :param pos_array: crystal to find forces of.
        :return: a vector of size 2N describing the x forces and y forces.
        """
        nion = self.nion
        x = pos_array[0:nion]
        y = pos_array[nion:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx ** 2 + dy ** 2)

        wr = self.geonium.to_theory_units('wall', self.geonium.trap.wall_f)
        wc = self.geonium.to_theory_units('cyclotron', self.geonium.cyc_f)
        m = self.geonium.to_theory_units('mass', self.geonium.m)
        q = self.geonium.q
        wall_order = self.geonium.trap.wall_order
        wc2 = (wall_order==2)*q*self.geonium.trap.wall_pot_expansion_in_rphi_coef(2)
        wc3 = (wall_order==3)*q*self.geonium.trap.wall_pot_expansion_in_rphi_coef(3)
        b = self.geonium.trap.b
        k_e = self.geonium._k_e

        # Calculate coulomb force on each ion
        with np.errstate(divide='ignore'):
            fc = np.where(rsep != 0., rsep ** (-2), 0)

        with np.errstate(divide='ignore', invalid='ignore'):
            fx = np.where(rsep != 0., np.float64((dx / rsep) * fc), 0)
            fy = np.where(rsep != 0., np.float64((dy / rsep) * fc), 0)

        # total force on each ion

        ftrapx = (-m*wr**2 - q*wc2 + q*b*wr + 2*wc) * x \
            - 4*q*wc3 * (x**3 + x*y**2) + 3*wc3*(x**2 - y**2)
        ftrapy = (-m*wr**2 - q*wc2 + q*b*wr - 2*wc2) * y \
            - 4*q*wc3 * (y**3 + y*x**2) - 6*wc3*x*y

        fx = -np.sum(fx, axis=1) + ftrapx
        fy = -np.sum(fy, axis=1) + ftrapy

        return np.array([fx, fy]).flatten()


class ModesTransverse:
    """
    Compute transverse modes.
    """
    def __init__(self, crystal: IonCrystal2d):
        """

        :param crystal: instance of class IonCrystal2d
        :return:
        """
        self.crystal = crystal
        pass

    def jij(self, odf_df):
        evec = 0
        eval = 0
        return ModesTransverse.__spin_spin_jij(odf_df, evec, eval)

    @staticmethod
    def __spin_spin_jij(odf_df, evec, eval):
        # mu is in rad/sec
        nions = len(eval)
        jij = np.zeros((nions, nions))
        for i in range(nions):
            for j in range(nions):
                jij[i,j] = np.sum(evec[:,i]*evec[:,j]/(odf_df**2-eval**2))
        return np.matrix(jij)

    @staticmethod
    def check_single_plane(eval):
        """Single-plane check

        :param eval: eigenfrequencies
        :return: bool
        """
        num_imag_eval = np.sum([np.imag(x)>0 for x in eval])
        if num_imag_eval > 0:
            print("found {} imaginary eigenvalues".format(num_imag_eval))
            return False
        return True

    def calc_axial_modes_simple(self, pos_array):
        """Solution for axial eigenvalues. Based on Adam Keith's 2011
        MATLAB code.

        :param pos_array: Position vector which defines the crystal to be analyzed.

        :return: Array of eigenvalues, Array of eigenvectors
        """
        m = self.m[0]
        N = int(pos_array.size / 2)
        A = np.empty((N, N))
        q = self.q
        k_e = self.k_e
        x = pos_array[0:N]
        y = pos_array[N:]
        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx ** 2 + dy ** 2)

        with np.errstate(divide='ignore'):
            rsep3 = np.where(rsep != 0., rsep**(-3), 0)

        A1 = np.diag((2 * q * self.Coeff[2] - k_e * q ** 2 * np.sum(rsep3, axis=0)))
        A2 = k_e * q ** 2 * rsep3
        A[:] = (A1 + A2) / m

        eval, evect = np.linalg.eig(A)
        eval = np.lib.scimath.sqrt(eval)
        return eval, evect

    def calc_axial_modes(self, pos_array):
        """
        Calculate the modes of axial vibration for a crystal defined
        by pos_array.

        THIS MAY NEED TO BE EDITED FOR NONHOMOGENOUS MASSES

        :param pos_array: Position vector which defines the crystal
                          to be analyzed.
        :return: Array of eigenvalues, Array of eigenvectors
        """

        x = pos_array[0:self.Nion]
        y = pos_array[self.Nion:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx ** 2 + dy ** 2)

        with np.errstate(divide='ignore'):
            rsep3 = np.where(rsep != 0., rsep ** (-3), 0)

        K = np.diag((-1 + 0.5 * np.sum(rsep3, axis=0)))
        K -= 0.5 * rsep3
        # Make first order system by making space twice as large
        Zn = np.zeros((self.Nion, self.Nion))
        eyeN = np.identity(self.Nion)
        Mmat = np.diag(self.md)
        Minv = np.linalg.inv(Mmat)
        firstOrder = np.bmat([[Zn, eyeN], [np.dot(Minv,K), Zn]])
        Eval, Evect = np.linalg.eig(firstOrder)

        # Convert 2N imaginary eigenvalues to N real eigenfrequencies
        ind = np.argsort(np.absolute(np.imag(Eval)))
        Eval = np.imag(Eval[ind])
        Eval = Eval[Eval >= 0]      # toss the negative eigenvalues
        Evect = Evect[:, ind]     # sort eigenvectors accordingly

        # Normalize by energy of mode
        for i in range(2*self.Nion):
            pos_part = Evect[:self.Nion, i]
            vel_part = Evect[self.Nion:, i]
            norm = vel_part.H*Mmat*vel_part - pos_part.H*K*pos_part

            with np.errstate(divide='ignore',invalid='ignore'):
                Evect[:, i] = np.where(np.sqrt(norm) != 0., Evect[:, i]/np.sqrt(norm), 0)
            #Evect[:, i] = Evect[:, i]/np.sqrt(norm)

        Evect = np.asarray(Evect)
        return Eval, Evect

    def is_plane_stable(self):
        """
        Checks to see if any of the axial eigenvalues in the current configuration of the crystal
        are equal to zero. If so, this indicates that the one-plane configuration is unstable
        and a 1-2 plane transistion is possible.

        :return: Boolean: True if no 1-2 plane transistion mode exists, false if it does
        (Answers: "is the plane stable?")

        """
        if self.hasrun is False:
            self.run()

        for x in self.axialEvals:
            if x == 0.0:
                return False

        return True

    def show_axial_evals(self, experimentalunits=False, flatlines=False):
        """
        Plots the axial eigenvalues vs mode number.
        :param experimentalunits:
        :return:
        """
        if self.axialEvals is []:
            print("Warning-- no axial eigenvalues found. Cannot show axial eigenvalues")
            return False

        if flatlines is True:
            fig = plt.figure(figsize=(8, 5))
            fig = plt.axes(frameon=True)
            # fig= plt.axes.get_yaxis().set_visible(False)
            fig.set_yticklabels([])

            if experimentalunits is False:
                fig = plt.xlabel("Eigenfrequency (Units of $\omega_z$)")
                fig = plt.xlim(min(self.axialEvals) * .99, max(self.axialEvals * 1.01))

                for x in self.axialEvals:
                    fig = plt.plot([x, x], [0, 1], color='black', )
            if experimentalunits is True:
                fig = plt.xlabel("Eigenfrequency (2 \pi* Hz)")
                fig = plt.xlim(min(self.axialEvalsE) * .99, max(self.axialEvalsE) * 1.01)
                for x in self.axialEvalsE:
                    fig = plt.plot([x, x], [0, 1], color='black')

            fig = plt.ylim([0, 1])
            # fig= plt.axes.yaxis.set_visible(False)
            fig = plt.title("Axial Eigenvalues for %d Ions, $f_{rot}=$%.1f kHz, and $V_{wall}$= %.1f V " %
                            (self.Nion, self.wrot / (2 * pi * 1e3), self.Cw * self.V0 / 1612))
            plt.show()
            return True

        fig = plt.figure()
        xvals = np.array(range(self.Nion))
        xvals += 1
        if experimentalunits is False:
            fig = plt.plot(xvals, sorted(self.axialEvals), "o")
            fig = plt.ylim((.97 * min(self.axialEvals), 1.01 * max(self.axialEvals)))
            fig = plt.ylabel("Eigenfrequency (Units of $\omega_z$)")
            fig = plt.plot([-1, self.Nion + 1], [1, 1], color="black", linestyle="--")

        else:
            fig = plt.plot(xvals, sorted(self.axialEvalsE), "o")
            fig = plt.ylim(.95 * min(self.axialEvalsE), 1.05 * max(self.axialEvalsE))
            fig = plt.ylabel("Eigenfrequency (Hz)")
            fig = plt.plot([-1, self.Nion + 1], [max(self.axialEvalsE), max(self.axialEvalsE)], color="black",
                           linestyle="--")

        fig = plt.xlabel("Mode Number")

        fig = plt.title("Axial Eigenvalues for %d Ions, $f_{rot}=$%.1f kHz, and $V_{wall}$= %.1f V " %
                        (self.Nion, self.wrot / (2 * pi * 1e3), self.Cw * self.V0 / 1612))
        fig = plt.xlim((0, self.Nion + 1))
        fig = plt.grid(True)
        fig = plt.show()
        return True

    def show_crystal_modes(self, pos_vect, Evects, modes):
        """
        For a given crystal, plots the crystal with colors based on the eigenvectors.

        :param pos_vect: the position vector of the current crystal
        :param Evects: the eigenvectors to be shown
        :param modes: the number of modes you wish to see
        """
        plt.figure(1)

        for i in range(modes):
            plt.subplot(modes, 1, i + 1, aspect='equal')
            plt.scatter(1e6 * pos_vect[0:self.Nion], 1e6 * pos_vect[self.Nion:], c=Evects[:, i], vmin=-.25, vmax=0.25,
                        cmap='RdGy')
            plt.xlabel('x position [um]')
            plt.ylabel('y position [um]')
            plt.axis([-200, 200, -200, 200])
        plt.tight_layout()

    def show_low_freq_mode(self):
        """
        Gets the lowest frequency modes and eigenvectors,
        then plots them, printing the lowest frequency mode.

        """
        num_modes = np.size(self.Evals)
        low_mode_freq = self.Evals[-1]
        low_mode_vect = self.Evects[-1]

        plt.scatter(1e6 * self.u[0:self.Nion], 1e6 * self.u[self.Nion:],
                    c=low_mode_vect, vmin=-.25, vmax=0.25, cmap='RdGy')
        plt.axes().set_aspect('equal')
        plt.xlabel('x position [um]', fontsize=12)
        plt.ylabel('y position [um]', fontsize=12)
        plt.axis([-300, 300, -300, 300])
        print(num_modes)
        print("Lowest frequency mode at {0:0.1f} kHz".format(float(np.real(low_mode_freq))))
        return 0

class ModesInPlane:
    """
    Compute in-plane modes for 2-d crystal
    """
    def __init__(self,
                 crystal: IonCrystal2d,
                 axial: ModesTransverse):
        """

        :param crystal: instance of class IonCrystal2d
        :param axial: instance of class ModesTransverse
        :return:
        """
        self.crystal = crystal

        self.axialEvals = []  # Axial eigenvalues
        self.axialEvects = []  # Axial eigenvectors
        self.planarEvals = []  # Planar eigenvalues
        self.planarEvects = []  # Planar Eigenvectors

        self.axialEvalsE = []  # Axial eigenvalues in experimental units
        self.planarEvalsE = []  # Planar eigenvalues in experimental units

        self.p0 = 0    # dimensionless potential energy of equilibrium crystal
        self.r = []
        self.rsep = []
        self.dx = []
        self.dy = []

    def calc_planar_modes(self, pos_array):
        """Calculate Planar Mode Eigenvalues and Eigenvectors

        THIS MAY NEED TO BE EDITED FOR NONHOMOGENOUS MASSES

        :param pos_array: Position vector which defines the crystal
                          to be analyzed.
        :return: Array of eigenvalues, Array of eigenvectors
        """

        V = -self.hessian_penning(pos_array)  # -Hessian
        Zn = np.zeros((self.Nion, self.Nion))
        Z2n = np.zeros((2 * self.Nion, 2 * self.Nion))
        offdiag = (2 * self.wr - self.wc) * np.identity(self.Nion)
        A = np.bmat([[Zn, offdiag], [-offdiag, Zn]])

        Mmat = np.diag(np.concatenate((self.md,self.md)))
        Minv = np.linalg.inv(Mmat)

        firstOrder = np.bmat([[Z2n, np.identity(2 * self.Nion)], [np.dot(Minv,V/2), A]])
        #mp.dps = 25
        #firstOrder = mp.matrix(firstOrder)
        #Eval, Evect = mp.eig(firstOrder)
        Eval, Evect = np.linalg.eig(firstOrder)
        # currently giving too many zero modes (increase numerical precision?)

        # make eigenvalues real.
        ind = np.argsort(np.absolute(np.imag(Eval)))
        Eval = np.imag(Eval[ind])
        Eval = Eval[Eval >= 0]      # toss the negative eigenvalues
        Evect = Evect[:, ind]    # sort eigenvectors accordingly

        # Normalize by energy of mode
        for i in range(4*self.Nion):
            pos_part = Evect[:2*self.Nion, i]
            vel_part = Evect[2*self.Nion:, i]
            norm = vel_part.H*Mmat*vel_part - pos_part.H*(V/2)*pos_part

            with np.errstate(divide='ignore'):
                Evect[:, i] = np.where(np.sqrt(norm) != 0., Evect[:, i]/np.sqrt(norm), 0)
            #Evect[:, i] = Evect[:, i]/np.sqrt(norm)

        # if there are extra zeros, chop them
        Eval = Eval[(Eval.size - 2 * self.Nion):]
        return Eval, Evect


class Visualize:
    def __init__(self,
                 crystal: IonCrystal2d,
                 axialm: ModesTransverse=None,
                 inplanem: ModesInPlane=None):
        """
        :param crystal: instance of class IonCrystal2d
        :param axialm: instance of class ModesTransverse
        :param inplanem: instance of class ModesInPlane
        :return:
        """
        self.crystal = crystal
        self.axialm = axialm
        self.inplanem = inplanem

    def show_crystal(self):
        """
        Makes a pretty plot of the crystal with a given position vector.

        :param pos_vect: The crystal position vector to be seen.
        """
        x, y = self.crystal.xy
        plt.plot(x, y, '.')
        plt.xlabel('x (um)')
        plt.ylabel('y (um)')
        plt.axes().set_aspect('equal')

        plt.show()

    @staticmethod
    def get_git_revision_short_hash():
        s = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
        return s.decode("utf-8").strip("\n")

    def visualize_ions(self, odf_df=100, to_disk=False):
        """Display summary of Mode Analysis Code results.

        :param ions: ModeAnalysis object
        :param odf_df: optical dipole force detuning (Hz)
        :param to_disk: write to disk (bool)
        :return: none
        """
        fig = plt.figure(0, figsize=(8,6), dpi=150)
        plt.clf()
        grid = gridspec.GridSpec(3,3)
        aximage = plt.subplot(grid[1,0])
        axhist = plt.subplot(grid[0,0:2])
        axtext = plt.subplot(grid[2,0])
        #axrt = plt.subplot(grid[:,2])
        axjij = plt.subplot(grid[1,1])
        # ax4 = plt.subplot2grid(gridt, (1, 0), colspan=2, rowspan=2)

        fig.suptitle("mode_analysis_code.py visualize_ions", fontsize=12)

        x, y = self.crystal.xy

        # ion image
        aximage.axvline(x=0, color="red")
        aximage.axhline(y=0, color="red")
        aximage.plot(x*1e6, y*1e6, '.')
        # ax4.set_xlabel('x position [um]')
        # ax4.set_ylabel('y position [um]')
        aximage.set_aspect('equal')
        aximage.axis([-300, 300, -300, 300])

        return

        ############################
        # calculate transverse modes
        ############################
        modes = ions.calc_axial_modes_simple(xy)
        crystal.check_single_plane(modes[0])
        f_mag = ions.wmag/2/np.pi
        f_com = modes[0][0]/2/np.pi
        f_modes = modes[0]/2/np.pi
        eval = modes[0]
        evec = np.transpose(modes[1])

        ################
        # mode histogram
        ################
        range=(f_mag/1e3, f_com/1e3)
        num_bins = 30
        hist, bins, patches = axhist.hist(f_modes/1e3, bins=num_bins )
        # ax1.set_xlabel("Mode Frequency (kHz)")
        # ax1.set_ylabel("Count in Bin")
        axhist.axvline(x=f_com/1e3, color="green")
        com_label = "{:.1f} kHz".format(f_com/1e3)
        axhist.text(x=f_com/1e3, y=np.max(hist)/2, s=com_label, rotation=90, size=12,
                        verticalalignment='center', horizontalalignment='right',
                        color = 'green', weight='bold')
        axhist.axvline(x=f_mag/1e3, color="green")
        mag_label = "{:.1f} kHz".format(f_mag/1e3)
        axhist.text(x=f_mag/1e3, y=np.max(hist)/2, s=mag_label, rotation=90, size=12,
                         verticalalignment='center', horizontalalignment='right',
                         color = 'green', weight='bold')

        ######################
        # eigenmodes
        ######################
        cnorm = evec[0]
        nhorz = 3
        nvert = 6
        subgridrt = gridspec.GridSpecFromSubplotSpec(nvert, nhorz,
                                                     subplot_spec=grid[:,2])
        ntoplot = min(nhorz*nvert, num_ions)
        for moden in np.arange(0, ntoplot):
            axt = plt.subplot(subgridrt[moden%nvert, moden//nvert])
            axt.scatter(x, y, c=evec[moden]/cnorm,
                         cmap=mpl.cm.bwr, linewidths=0)
            xyr = np.max(x)*1.25
            axt.set_xlim(-xyr, xyr)
            axt.set_ylim(-xyr, xyr)
            axt.axis('off')
            axt.set_aspect('equal')
            slabel = "{}".format(moden)
            axt.text(x=0, y=xyr, s=slabel, rotation=0, size=12,
                         verticalalignment='center', horizontalalignment='center',
                         color = 'black', weight='bold')
            slabel = "{:.1f} kHz".format(eval[moden]/2/np.pi/1e3)
            axt.text(x=0, y=-xyr*1.1, s=slabel, rotation=0, size=6,
                         verticalalignment='center', horizontalalignment='center',
                         color = 'black', weight='bold')

        ###########################
        # Jij as Figure 2 of Nature
        ###########################
        arb_norm = 1e-13
        dij, dij_row = calc_dij(x, y)
        w_odf = (f_com + odf_df)*2*np.pi
        Jij = spinspinJ(w_odf, evec, eval)/arb_norm  # Jij is np.matrix
        # don't count ions' self-interaction
        np.fill_diagonal(Jij, 0)  # in-place modification of J
        Jij_row = np.squeeze(np.array(np.reshape(Jij, (num_ions**2,1))))

        # subsample for plotting [start:stop:step]
        axjij.plot(dij_row*1e6, Jij_row, '+')
        axjij.set_yscale('log')
        axjij.set_xscale('log')
        # ax3.set_xlabel("dij (um)")
        # ax3.set_ylabel("Jij (a.u.)")
        axjij.set_xlim(10, 1000)
        axjij.set_ylim(np.min(Jij)/5,1)
        axjij.grid(b=True)

        #################
        # textual summary
        #################
        s = ""
        s += "S={}   N={}\n".format(ions.shells, ions.Nion)
        s += "Vtrap={} V [end, mid, center]\n".format(ions.Vtrap)
        s += "fwall={:.1f} kHz  \n"\
            .format(ions.wrot/2/np.pi/1e3)
        s += "\n"
        s += "f_com={:.1f} kHz  f_mag={:.1f} kHz\n".format(f_com/1e3, f_mag/1e3)
        s += "f_cyc={:.1f} MHz  mode_bw={:.1f} kHz\n"\
            .format(ions.wcyc/2/np.pi/1e6, (f_com-np.min(f_modes))/1e3)
        s += "min(dij)={:.1f} um   diameter={:.1f} um\n"\
            .format(calc_dij_min(x, y)*1e6, np.max(dij_row)*1e6)
        s += "ODF_detuning={:.1f} kHz\n".format(odf_df/1e3)
        s += "\n"
        s += "github.com/nistpenning/calc : {}"\
            .format(get_git_revision_short_hash())
        axtext.text(0, 0, s, size=8)
        axtext.set_axis_off()
        axtext.set_axis_bgcolor('white')
        plt.show()

        ###############
        # write to disk
        ###############
        if to_disk:
            tstamp = int(time.time())

            # write image
            fname = "mode_code_{}.png".format(tstamp)
            plt.savefig(fname, bbox_inches='tight')

            # write configuration data
            fname = "mode_code_config_{}.csv".format(tstamp)
            fh = open(fname,'w')
            fh.write(s)
            fh.close()

            # write data
            fmts='%1.4e'
            fname = "mode_code_{}_xy.csv".format(tstamp)
            np.savetxt(fname, np.array([x, y]), delimiter=", ",
                       header=s, fmt=fmts)
            fname = "mode_code_{}_eigen_val.csv".format(tstamp)
            np.savetxt(fname, eval/2/np.pi, delimiter=", ",
                       header=s, fmt=fmts)
            fname = "mode_code_{}_eigen_vect.csv".format(tstamp)
            np.savetxt(fname, evec, delimiter=", ",
                       header=s, fmt=fmts)
            fname = "mode_code_{}_Jij.csv".format(tstamp)
            np.savetxt(fname, Jij, delimiter=", ",
                       header=s, fmt=fmts)

