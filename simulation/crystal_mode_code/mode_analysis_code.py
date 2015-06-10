from __future__ import division, with_statement
from scipy.constants import pi
import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

__author__ = 'sbt'

# -*- coding: utf-8 -*-
"""
Contains the ModeAnalysis class, which can simulate the positions of ions in a crystal
of desired size. The class contains methods for the generation of a crystal,
relaxation to a minimum potential energy state, and determination of axial and (eventually) planar modes of motion
by methods derived by Wang, Keith, and Freericks in 2013.
"""


class ModeAnalysis:
    def __init__(self, shells=4, Vtrap=[0.0, -1750.0, -2000.0], Ctrap=1.0,
                 fz=1000, B=4.4588, frot=60., Vwall=1., wall_order=2, mult=1e14, quiet=True):
        """

        Initialization parameters:
        shells:     integer, number of shells to instantiate the plasma with
        Vrap:       array of 3 elements, defines the [end, middle, center] voltages on the trap electrodes.
        Ctrap:      float, constant coefficient on trap potentials
        B:          float, defines strength of axial magnetic field.
        frot:       float, frequency of rotation
        Vwall:      float, strength of wall potential in volts
        wall_order: integer, defines the order of the rotating wall potential
        mult:       float, mutliplicative factor for simplifying numerical calculations
        """
        self.Evects = []  # Fill list of eigenvectors
        self.Evals = []  # Fill list of eigenfrequencies
        self.dens = 0
        self.avg = 0
        self.quiet = quiet
        # Initialize basic variables such as physical constants
        self.shells = shells
        self.Nion = 1 + 6 * np.sum(range(1, shells + 1))
        self.u0 = np.empty(2 * self.Nion)  # for array of ion positions first half is x, last is y
        self.u = np.empty(2 * self.Nion)  # for array of ion positions first half is x, last is y
        # self.scale = 0
        self.q = 1.602E-19
        self.m_Be = 1.5E-26
        self.k_e = 8.9875517873681764E9
        # trap definitions
        self.B = B
        self.m = self.m_Be * np.ones(self.Nion)
        self.wcyc = self.q * B / self.m_Be  # Cyclotron frequency
        self.C = Ctrap * np.array([[0.0756, 0.5157, 0.4087],
                                   [-0.0001, -0.005, 0.005],
                                   [1.9197e3, 3.7467e3, -5.6663e3],
                                   [0.6738e7, -5.3148e7, 4.641e7]])  # axial trap coefficients; see Teale's final paper
        self.relec = 0.01  # rotating wall electrode distance in meters
        self.Vtrap = np.array(Vtrap)  # [Vend, Vmid, Vcenter] for trap electrodes
        self.Coeff = np.dot(self.C, self.Vtrap)  # Determine the 0th, first, second, and fourth order
        #                                          potentials at trap center

        self.wz = np.sqrt(2 * self.q * self.Coeff[2] / self.m_Be)  # Compute axial frequency
        self.wrot = 2 * pi * frot * 1e3  # Rotation frequency in units of angular frequency
        self.wmag = 0.5 * (self.wcyc - np.sqrt(self.wcyc ** 2 - 2 * self.wz ** 2))

        self.V0 = (0.5 * self.m_Be * self.wz ** 2) / self.q  # Find quadratic voltage at trap center
        self.Vw = self.V0 * 0.045 * Vwall / 1000  # V_T divides in MATLAB
        self.mult = mult  # store multiplicative factor

        # wall order
        if wall_order == 2:
            self.Cw2 = self.q * Vwall * 1612
            self.Cw3 = 0
        if wall_order == 3:
            self.Cw2 = 0
            self.Cw3 = self.q * Vwall * 3e4

    def run(self):
        """
        Generates a crystal by the find_scalled_lattice_guess method,
        adjusts it into an eqilibirium position by find_eq_pos method,
        and then computes the eigenvalues and eigenvectors of the axial modes by calc_axial_modes.

        Sorts the eigenvalues and eigenvectors and stores them in self.Evals, self.Evects.
        Stores the radial separations as well.
        """
        mins = 10e-6
        res = 0.1e-6
        if self.wmag > self.wrot:
            print("Warning: Rotation frequency below magnetron frequency of {0:.1f}".format(float(self.wmag / 2 * pi)))
            return 0
        self.u0[:] = self.find_scaled_lattice_guess(mins, res)
        self.u = self.find_eq_pos(self.u0)

        self.Evals, self.Evects = self.calc_axial_modes(self.u)

        # sort arrays
        sort_ind = np.argsort(np.abs(self.Evals))[::-1]
        self.Evals = self.Evals[sort_ind] / (2 * pi * 1e3)  # units of kHz
        self.Evects = self.Evects[:, sort_ind]
        self.r, dx, dy, self.rsep = self.find_radial_separation(self.u)

    def find_scaled_lattice_guess(self, mins, res):
        """
        Will generate a 2d hexagonal lattice based on the shells intialiization parameter.
        Guesses initial minimum separation of mins and then increases spacing until a local minimum of
        potential energy is found.
        """

        # Make a 2d lattice; u represents the position
        uthen = self.generate_2D_hex_lattice(mins)
        # Figure out the lattice's initial potential energy
        pthen = self.pot_energy(uthen)

        # Iterate through the range of minimum spacing in steps of res/resolution
        for scale in np.arange(mins + res, (mins / res) * mins, res):
            # Quickly make a 2d hex lattice; perhaps with some stochastic procedure?
            uguess = self.generate_2D_hex_lattice(scale)
            # Figure out the potential energy of that newly generated lattice
            pnow = self.pot_energy(uguess)

            # And if the program got a lattice that was less favorably distributed, conclude
            # that we had a pretty good guess and return the lattice.
            if pnow >= pthen:
                # print "find_scaled_lattice: Minimum found"
                # print "initial scale guess: " + str(scale)
                # self.scale = scale
                return uthen
            # If not, then we got a better guess, so store the energy score and current arrangement
            # and try again for as long as we have mins and resolution to iterate through.
            uthen = uguess
            pthen = pnow
        # If you're this far it means we've given up
        # self.scale = scale
        # print "find_scaled_lattice: no minimum found, returning last guess"
        return uthen

    # Make a 2d lattice
    def generate_2D_hex_lattice(self, scale=1):
        # Start out with the initial
        posvect = np.array([0.0, 0.0])  # always a point [0,0]

        for sh in range(1, self.shells + 1):
            posvect = np.append(posvect, self.add_hex_shell(sh))  # Inefficiency; instead of appending could fill out
        posvect *= scale
        # Posvect is refolded
        posvect = np.transpose(np.reshape(posvect, (self.Nion, 2)))
        # Then refolded again? and returned
        return np.reshape(posvect, 2 * self.Nion)

    # A slave function used to append shells onto a position vector
    def add_hex_shell(self, s):
        a = list(range(s, -s - 1, -1))  # Change 1
        a.extend(-s * np.ones(s - 1))
        a.extend(range(-s, s + 1))
        a.extend(s * np.ones(s - 1))

        b = list(range(0, s + 1))
        b.extend(s * np.ones(s - 1))
        b.extend(range(s, -s - 1, -1))
        b.extend(-s * np.ones(s - 1))
        b.extend(range(-s, 0))

        x = np.sqrt(3) / 2.0 * np.array(b)
        y = 0.5 * np.array(b) + np.array(a)
        pair = np.column_stack((x, y)).flatten()

        return pair

    # Get the potential energy from a positional array
    def pot_energy(self, pos_array):
        # Frequency of rotation, mass and the number of ions in the array
        wr = self.wrot
        m = self.m[0]
        q = self.q
        B = self.B
        k_e = self.k_e
        N = int(pos_array.size / 2)

        # the x positions are the first N elements of the position array
        x = pos_array[0:N]
        # The y positions are the last N elements of the position array
        y = pos_array[N:]

        # dx flattens the array into a row and 'normalizes' by subtracting itself to get some zeroes.
        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        # rsep is the distances between
        rsep = np.sqrt(dx ** 2 + dy ** 2)

        # dxsq = np.array([(i-x)**2 for i in x])
        # dysq = np.array([(i-y)**2 for i in y])

        # dxsq = [(i-x)**2 for i in x]
        # dysq = [(i-y)**2 for i in y]

        # rsep = np.sqrt(dxsq + dysq)

        with np.errstate(divide='ignore'):
            Vc = np.where(rsep != 0., 1/rsep, 0)

        # One half times the rotational force, the charge times the coeff,
        V = 0.5 * (-m * wr ** 2 - q * self.Coeff[2] + q * B * wr) * np.sum((x ** 2 + y ** 2)) \
            - q * self.Coeff[3] * np.sum((x ** 2 + y ** 2) ** 2) \
            + np.sum(self.Cw2 * (x ** 2 - y ** 2)) \
            + np.sum(self.Cw3 * (x ** 3 - 3 * x * y ** 2)) \
            + 0.5 * k_e * q ** 2 * np.sum(Vc)
        return self.mult * V

    def force_penning(self, pos_array):
        wr = self.wrot
        m = self.m[0]
        q = self.q
        B = self.B
        N = int(pos_array.size / 2)

        x = pos_array[0:N]
        y = pos_array[N:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx ** 2 + dy ** 2)

        # Calculate coulomb force on each ion
        with np.errstate(divide='ignore'):
            Fc = np.where(rsep != 0., rsep**(-2), 0)

        with np.errstate(divide='ignore', invalid='ignore'):
            fx = np.where(rsep != 0., np.float64((dx / rsep) * Fc), 0)
            fy = np.where(rsep != 0., np.float64((dy / rsep) * Fc), 0)

        # total force on each ion

        Ftrapx = (-m * wr ** 2 - q * self.Coeff[2] + q * B * wr + 2 * self.Cw2) * x \
            - 4 * q * self.Coeff[3] * (x ** 3 + x * y ** 2) + 3 * self.Cw3 * (x ** 2 - y ** 2)
        Ftrapy = (-m * wr ** 2 - q * self.Coeff[2] + q * B * wr - 2 * self.Cw2) * y \
            - 4 * q * self.Coeff[3] * (y ** 3 + y * x ** 2) - 6 * self.Cw3 * x * y

        # Ftrap =  (m*w**2 + q*self.V0 - 2*q*self.Vw - q*self.B* w) * pos_array
        Fx = -(self.k_e * self.q ** 2) * np.sum(fx, axis=1) + Ftrapx
        Fy = -(self.k_e * self.q ** 2) * np.sum(fy, axis=1) + Ftrapy

        Fx[np.abs(Fx) < 1e-24] = 0
        Fy[np.abs(Fy) < 1e-24] = 0

        Fx *= self.mult
        Fy *= self.mult
        return np.array([Fx, Fy]).flatten()

    def find_eq_pos(self, u0):
        fun_tolerance = 1e-37

        out = optimize.minimize(self.pot_energy, u0, method='BFGS', jac=self.force_penning,
                                options={'gtol': fun_tolerance, 'disp': not self.quiet})
        return out.x

    def calc_axial_modes(self, pos_array):
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

        Eval, Evect = np.linalg.eig(A)

        Eval = np.lib.scimath.sqrt(Eval)

        return Eval, Evect

    def show_crystal(self, pos_vect):
        plt.plot(1e6 * pos_vect[0:self.Nion], 1e6 * pos_vect[self.Nion:], '.')
        plt.xlabel('x position [um]')
        plt.ylabel('y position [um]')
        plt.axes().set_aspect('equal')
        plt.axis([-300, 300, -300, 300])

        plt.show()

    def show_crystal_modes(self, pos_vect, Evects, modes):
        plt.figure(1)

        for i in range(modes):
            plt.subplot(modes, 1, i + 1, aspect='equal')
            plt.scatter(1e6 * pos_vect[0:self.Nion], 1e6 * pos_vect[self.Nion:], c=Evects[:, i], vmin=-.25, vmax=0.25,
                        cmap='RdGy')
            plt.xlabel('x position [um]')
            plt.ylabel('y position [um]')
            plt.axis([-200, 200, -200, 200])
        plt.tight_layout()

    def get_low_freq_mode(self):
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

    @staticmethod
    def nan_to_zero(my_array):
        my_array[np.isinf(my_array) | np.isnan(my_array)] = 0
        return my_array

    @staticmethod
    def save_positions(u):
        np.savetxt("py_u.csv", u, delimiter=",")

    @staticmethod
    def crystal_spacing_fit(r, offset, curvature):
        return np.sqrt(2 / (np.sqrt(3) * offset * np.sqrt(1 - (r * curvature) ** 2)))

    @staticmethod
    def find_radial_separation(pos_array):
        N = int(pos_array.size / 2)
        x = pos_array[0:N]
        y = pos_array[N:]
        r = np.sqrt(x ** 2 + y ** 2)

        sort_ind = np.argsort(r)

        r = r[sort_ind]
        x = x[sort_ind]
        y = y[sort_ind]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx ** 2 + dy ** 2)

        return r, dx, dy, rsep


########################################################################################

if __name__ == "__main__":
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculationConsistency)
    # unittest.TextTestRunner(verbosity=1).run(suite)
    shellcounts = [4, 5]
    transistionfreq = []
    ions = []
    for s in shellcounts:
        for w in np.linspace(185, 220, 10):
            a = ModeAnalysis(shells=s, Vtrap=[0.0, -1750.0, -2000.0], Ctrap=1.0, frot=w, Vwall=2.13, wall_order=2)
            a.run()
            print(a.Nion)
            Evals = a.Evals
            # print(Evals)
            print("Trying w=", w)
            # print("Calculated axial freq:",a.wz/(2*pi))
            if type(Evals) is type(1):
                print("Integer eval obtained:", Evals)
                # else:
                # print(Evals[0])
            if [x for x in Evals if np.imag(x) != 0] != []:
                ions.append(a.Nion)
                transistionfreq.append(w)
                # print([x for x in Evals if np.imag(x) != 0])
                break
            print("------------")
    print(transistionfreq)
    a = plt.figure()
    a = plt.plot(shellcounts, transistionfreq)
    plt.show()
    # #"""
    # #print("Checkpoint Alpha")
    # #plt.close()
    # #plt.plot((Evals,Evals),(np.ones(np.size(Evals)),np.zeros(np.size(Evals))),linestyle="-",color='black')
    # #plt.axis([800,1.550,0,1.1])
    # #plt.ylabel('arb')
    # #plt.xlabel('Mode Frequency [kHz]')
    # #a.show_crystal_modes(a.u, a.Evects, 3)
    # #print(np.max(Evals))
    # #print(np.min(Evals))
    # #plt.close()
    # a.show_crystal(a.u)
    # # print(a.Nion)
    #
    # #   get the data
    # r = np.zeros(np.shape(a.r))
    # rsep = np.zeros(np.shape(a.rsep))
    # r[:] = a.r
    # rsep[:] = np.transpose(a.rsep) / a.scale
    # #    print("Checkpoint Charlie")
    # a.get_low_freq_mode()
    # #   process into a single 1D array -- use braodcasting to make a r into repeating
    # #   2D array (faster than a loop)
    # r_full = r + np.zeros(np.shape(rsep))
    # rsep_flat = rsep.flatten()
    # r_flat = r_full.flatten()
    #
    # #   create mask for points
    # rsep_flat[rsep_flat > 1.4] = np.nan
    # rsep_flat[rsep_flat < 0.1] = np.nan
    # mask = np.isnan(rsep_flat)
    # r_fit = r_flat[~mask]
    # rsep_fit = rsep_flat[~mask]
    #
    # popt, pcov = optimize.curve_fit(a.crystal_spacing_fit, r_fit, rsep_fit,
    #                                 p0=[1.25, 2000])
    #
    # pred = a.crystal_spacing_fit(r, *popt)
    # plt.figure()
    # plt.plot(r_flat / a.scale, rsep_flat, 'o')
    # plt.plot(r / a.scale, pred, 'r')
    # plt.ylim(0, 2)
    # plt.show()
    # print(popt)
