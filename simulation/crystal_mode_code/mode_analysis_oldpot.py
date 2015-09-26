# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 14:12:58 2015

@author: ack
Questions:
    Still stuck on new potentials
    What is Vwall vs Vw?

Future Changes:
    Allow for any number of ions
"""


from scipy.constants import pi
import numpy as np
import scipy.constants as cons

import scipy.optimize as optimize
import matplotlib.pyplot as plt

class ModeAnalysis():
    """Zero temperature normal mode calculator"""
    
    def __init__(self, shells=4, Vtrap = [0.0,-1750.0,-2000.0], Ctrap = 1.0, 
                 fz=1000, B=4.4588, frot=60, Vwall=1, wall_order=2, mult=1):
        """Initialize Trap Parameters"""
        
        # Physical constants
        q = cons.elementary_charge        
        epsilon0 = cons.epsilon_0
        k_e = 1/(4*pi*epsilon0) #Coulombs constant
        amu = 1.66057e-27
        self.m_Be = 9.012182*amu
        
        self.shells = shells
        self.Nion = 1 + 6* np.sum(np.arange(shells+1)) # number of ions
        
        # ion position arrays first half is x, last is y
        self.u0 = np.empty(2*self.Nion) 
        self.u = np.empty(2*self.Nion) 
        self.scale = 0
        
        # Trap parameters
        self.C = Ctrap * np.array([[0.0756,0.5157,0.4087],
                   [-0.0001,-0.005,0.005],
                   [1.9197e3, 3.7467e3, -5.6663e3],
                   [0.6738e7, -5.3148e7, 4.641e7]]) #axial trap coefficients
        self.relec = 0.01  # rot wall electrode distance in meters
        self.Vtrap = np.array(Vtrap)  # Vend, Vmid, Vcenter
        self.Coeff = np.dot(self.C, self.Vtrap)

        self.B = B
        self.m = self.m_Be * np.ones(self.Nion)
        self.wcyc = q*self.B/self.m_Be
        #self.wz = np.sqrt(2*q*self.Coeff[2]/self.m_Be)  # why?
        self.wz = 4.9951e6  # check old stuff
        self.wmag = 0.5*(self.wcyc-np.sqrt(self.wcyc**2-2*self.wz**2))
        self.wrot = 2*pi*frot*1e3

        self.V0 = (0.5*self.m_Be*self.wz**2)/q
        self.Vw = self.V0*0.045*Vwall/1000  # V_T divides in MATLAB
        self.mult = mult

        # Wall order
        if wall_order == 2:
            self.Cw2 = q*Vwall*1612
            self.Cw3 = 0
        if wall_order == 3:
            self.Cw2 = 0
            self.Cw3 = q*Vwall*3e4

        # Convert to dimensionless quantities
        self.m = self.m_Be * np.ones(self.Nion)    # dimensionless mass
        self.l0 = ((k_e*q**2)/(q*self.V0))**(1/3)  # dimensionless length
        self.t0 = 1/self.wz                        # dimensionless time
        self.v0 = self.l0/self.t0                  # dimensionless velocity
        self.wr = self.wrot/self.wz                # dimensionless rotation
        self.wc = self.wcyc/self.wz                # dimensionless cyclotron

    def run(self, quiet=False):
        """Calculate Everything"""

        if self.wmag > self.wrot:
            print("Warning: Rotation frequency below magnetron frequency of {0.1f}".format(self.wmag))
        #self.u0[:] = self.find_scaled_lattice_guess(mins, res)
        self.u0 = generate_2D_hex_lattice(self.shells)            
        #self.show_crystal(self.u0)
        #print("Calculating Equilibrium Positions")
        self.u = self.find_eq_pos(self.u0)
        
        #self.show_crystal(self.u)
        #print "Calculate transverse axial modes"
        self.axialEvals, self.axialEvects = self.calc_axial_modes(self.u)
        self.planarEvals, self.planarEvects = self.calc_planar_modes(self.u)
        
        #sort arrays
        #sort_ind = np.argsort(np.abs(self.Evals))[::-1]
        #self.Evals = self.Evals[sort_ind]/(2*pi*1e3) #units of kHz
        #self.Evects = self.Evects[:,sort_ind]
        self.r,dx,dy,self.rsep = self.find_radial_separation(self.u)
        
        #plt.plot((self.Evals,self.Evals),(np.ones(np.size(self.Evals)),np.zeros(np.size(self.Evals))),linestyle="-",color='black')
        #plt.axis([800,1010,0,1.1])
        
        #self.show_crystal_modes(self.u, self.Evects,3)
        
    def pot_energy(self, pos_array):
        """Calculate Potential Energy of crystal"""

        x = pos_array[0:self.Nion]
        y = pos_array[self.Nion:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx**2 + dy**2)

        Vc = 1/rsep
        Vc[np.isinf(Vc) | np.isnan(Vc)] = 0

        V = -np.sum((self.m*self.wr**2 - self.wr*self.wc + 0.5)*(x**2 + y**2))\
            + (self.Vw/self.V0)*np.sum(x**2 + y**2) + 0.5*np.sum(Vc)

        return self.mult*V

    def force_penning(self, pos_array):
        """Calculate (positive) gradient of potential"""

        x = pos_array[0:self.Nion]
        y = pos_array[self.Nion:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx**2 + dy**2)

        # Calculate coulomb force on each ion
        Fc = rsep**(-2)
        Fc[np.isinf(Fc) | np.isnan(Fc)] = 0
        fx = np.float64((dx/rsep)*Fc)
        fx[np.isinf(fx) | np.isnan(fx)] = 0
        fy = np.float64((dy/rsep)*Fc)
        fy[np.isinf(fy) | np.isnan(fy)] = 0

        # Total force on each ion
        Ftrapx = -2*(self.m*self.wr**2 - self.wr*self.wc + 0.5 -
                     self.Vw/self.V0)*x
        Ftrapy = -2*(self.m*self.wr**2 - self.wr*self.wc + 0.5 +
                     self.Vw/self.V0)*y

        Fx = -np.sum(fx, axis=1) + Ftrapx
        Fy = -np.sum(fy, axis=1) + Ftrapy

        Fx = self.mult*Fx
        Fy = self.mult*Fy
        #print(np.concatenate((Fx, Fy)))
        #return np.concatenate((Fx, Fy))
        return np.array([Fx,Fy]).flatten()

    def hessian_penning(self, pos_array):
        """Calculate Hessian of potential"""

        H = np.empty((2*self.Nion, 2*self.Nion))

        x = pos_array[0:self.Nion]
        y = pos_array[self.Nion:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx**2 + dy**2)
        rsep5 = self.nan_to_zero(rsep**-5)

        dxsq = dx**2
        dysq = dy**2

        # X derivatives, Y derivatives for alpha != beta
        Hxx = np.mat((rsep**2 - 3*dxsq)*rsep5)
        Hyy = np.mat((rsep**2 - 3*dysq)*rsep5)

        # Above, for alpha == beta
        Hxx += np.mat(np.diag(-2*(self.m*self.wr**2 - self.wr*self.wc + 1/2 -
                              self.Vw/self.V0) -
                              np.sum((rsep**2 - 3*dxsq)*rsep5, axis=0)))
        Hyy += np.mat(np.diag(-2*(self.m*self.wr**2 - self.wr*self.wc + 1/2 +
                              self.Vw/self.V0) -
                              np.sum((rsep**2 - 3*dxsq)*rsep5, axis=0)))

        # Mixed derivatives
        Hxy = np.mat(-3*dx*dy*rsep5)
        Hxy += np.mat(np.diag(3*np.sum(dx*dy*rsep5, axis=0)))

        H = np.bmat([[Hxx, Hxy], [Hxy, Hyy]])
        H = np.asarray(H)
        return self.mult*H

    def find_eq_pos(self, u0):
        fun_tolerance = 1e-37
        
        #return optimize.fmin_bfgs(self.pot_energy, u0, maxiter = 1000)
        
        #out = optimize.fmin_ncg(self.pot_energy, u0, self.force_penning, fhess=self.hessian_penning, 
        #                        maxiter=1000, ftol=fun_tolerance)
        
        #guessadd = np.array([4e-6*np.random.random() for i in u0])
                                
#        out = optimize.minimize(self.pot_energy, u0, method='Newton-CG', jac=self.force_penning, 
#                                hess=self.hessian_penning, tol=fun_tolerance, options={'xtol': 1e-37,'disp': True})
        self.quiet = False
        
        out = optimize.minimize(self.pot_energy, u0, method='Newton-CG', jac=self.force_penning, hess=self.hessian_penning,
                                options={'gtol': fun_tolerance, 'disp': not self.quiet})
        #out = optimize.minimize(self.pot_energy, u0, method='BFGS', jac=self.force_penning,
        #                        options={'gtol': fun_tolerance, 'disp': not self.quiet})
        
        return out.x

    def calc_axial_modes(self, pos_array):
        """Calculate Axial Mode Eigenvalues and Eigenvectors

        Assumes all ions same mass"""
        # A = np.empty((self.n, self.n))

        x = pos_array[0:self.Nion]
        y = pos_array[self.Nion:]

        dx = x.reshape((x.size, 1)) - x
        dy = y.reshape((y.size, 1)) - y
        rsep = np.sqrt(dx**2 + dy**2)
        rsep3 = self.nan_to_zero(rsep**-3)

        A = np.diag((1 - 0.5*np.sum(rsep3, axis=0)))
        A += 0.5*rsep3

        Eval, Evect = np.linalg.eig(A)
        Eval = np.lib.scimath.sqrt(Eval)
        return Eval, Evect

    def calc_planar_modes(self, pos_array):
        """Calculate Planar Mode Eigenvalues and Eigenvectors

        Assumes all ions same mass"""

        V = -self.hessian_penning(pos_array)  # -Hessian
        Zn = np.zeros((self.Nion, self.Nion))
        Z2n = np.zeros((2*self.Nion, 2*self.Nion)) 
        offdiag = (2*self.wr-self.wc)*np.identity(self.Nion)
        A = np.bmat([[Zn, offdiag], [-offdiag, Zn]])

        firstOrder = np.bmat([[Z2n, np.identity(2*self.Nion)], [V/2, A]])
        #firstOrder = mp.matrix(firstOrder)
        #Eval, Evect = mp.eig(firstOrder)

        Eval, Evect = np.linalg.eig(firstOrder)
        #Eval = np.lib.scimath.sqrt(Eval)
        return Eval, Evect
     
    
        
    def find_radial_separation(self, pos_array):
        N = int(pos_array.size/2)
        x = pos_array[0:N]
        y = pos_array[N:]
        r = np.sqrt(x**2+y**2)
        
        sort_ind = np.argsort(r)
        r = r[sort_ind]
        x = x[sort_ind]
        y = y[sort_ind]      
        
        dx = x.reshape((x.size,1)) - x
        dy = y.reshape((y.size,1)) - y
        rsep = np.sqrt(dx**2 + dy**2)
        
        return r,dx,dy,rsep
        
    def nan_to_zero(self, my_array):
        my_array[np.isinf(my_array) | np.isnan(my_array)] = 0
        return my_array
        
    def show_crystal(self, pos_vect):
        plt.plot(pos_vect[0:self.Nion],pos_vect[self.Nion:],'.')
        #ax.set_aspect('equal')
        plt.xlabel('x position [um]')
        plt.ylabel('y position [um]')
        plt.axes().set_aspect('equal')
        #plt.axis([-5,5,-5,5])     
        
        plt.show()
        
    def show_crystal_modes(self, pos_vect, Evects, modes):
        plt.figure(1)
        
        for i in range(modes):
            plt.subplot(modes,1,i+1,aspect='equal')
            plt.scatter(1e6*pos_vect[0:self.Nion],1e6*pos_vect[self.Nion:],c=Evects[:,i], vmin=-.25, vmax=0.25, cmap='RdGy')
            plt.xlabel('x position [um]')
            plt.ylabel('y position [um]')
            plt.axis([-200,200,-200,200])            
        plt.tight_layout()
    
    def get_low_freq_mode(self):
        num_modes = np.size(self.Evals)
        low_mode_freq = self.Evals[-1]
        low_mode_vect = self.Evects[-1]
        
        plt.scatter(1e6*self.u[0:self.Nion],1e6*self.u[self.Nion:],
                    c=low_mode_vect, vmin=-.25, vmax=0.25, cmap='RdGy')
        plt.axes().set_aspect('equal')
        plt.xlabel('x position [um]', fontsize=12)
        plt.ylabel('y position [um]', fontsize=12)
        plt.axis([-300,300,-300,300])
        print(num_modes)
        print("Lowest frequency mode at {0:0.1f} kHz".format(float(np.real(low_mode_freq))))
        return 0

    def save_positions(self, u):
        np.savetxt("py_u.csv", u, delimiter = ",")

    def crystal_spacing_fit(self, r, offset, curvature):
        return np.sqrt(2/(np.sqrt(3)*offset*np.sqrt(1-(r*curvature)**(2))))


def generate_2D_hex_lattice(shells=1, scale=1):
    """Generate closed shell hexagonal lattice"""
    posvect = np.array([0.0, 0.0])  # center ion at [0,0]

    for s in range(1, shells+1):
        posvect = np.append(posvect, add_hex_shell(s))
    posvect = scale * posvect
    return np.hstack((posvect[0:-1:2], posvect[1:-1:2], posvect[-1]))


def add_hex_shell(s):
    """Generate position for full shell s"""
    a = list(range(s, -s-1, -1))
    a.extend(-s*np.ones(s-1))
    a.extend(list(range(-s, s+1)))
    a.extend(s*np.ones(s-1))

    b = list(range(0, s+1))
    b.extend(s*np.ones(s-1))
    b.extend(list(range(s, -s-1, -1)))
    b.extend(-s*np.ones(s-1))
    b.extend(list(range(-s, 0)))

    x = np.sqrt(3)/2.0 * np.array(b)
    y = 0.5 * np.array(b) + np.array(a)
    pair = np.column_stack((x, y)).flatten()
    return pair

                
if __name__ == "__main__":
    a = ModeAnalysis(shells=8 ,Vtrap=[-0.0,-1750.0,-2000.0], Ctrap = 1.0, frot=180.0, Vwall= 200, wall_order=2)
    a.run()
    Evals = a.axialEvals
    Evect = a.axialEvects
    a.show_crystal(a.u)
#    plt.close()
#    plt.plot((Evals,Evals),(np.ones(np.size(Evals)),np.zeros(np.size(Evals))),linestyle="-",color='black')
#    plt.axis([800,1.550,0,1.1])
#    plt.ylabel('arb')
#    plt.xlabel('Mode Frequency [kHz]')
#    a.show_crystal_modes(a.u, a.Evects, 3)
#    print(np.max(Evals))
#    print(np.min(Evals))
#    plt.close()
#    a.show_crystal(a.u)
##    print(a.n)
#    
##   get the data    
#    r = np.zeros(np.shape(a.r))
#    rsep = np.zeros(np.shape(a.rsep))
#    r[:] = a.r
#    rsep[:]= np.transpose(a.rsep)/a.scale
#
##   process into a single 1D array -- use braodcasting to make a r into repeating 
##   2D array (faster than a loop)
#    r_full = r + np.zeros(np.shape(rsep)) 
#    rsep_flat = rsep.flatten()
#    r_flat = r_full.flatten()
#    
##   create mask for points
#    rsep_flat[rsep_flat>1.4] = np.nan
#    rsep_flat[rsep_flat<0.1] = np.nan
#    mask = np.isnan(rsep_flat)
#    r_fit = r_flat[~mask]
#    rsep_fit = rsep_flat[~mask]
#    
#    popt,pcov = optimize.curve_fit(a.crystal_spacing_fit, r_fit, rsep_fit,
#                                    p0=[1.25,2000])
#    
#    pred = a.crystal_spacing_fit(r, *popt)
#    
#    plt.plot(r_flat/a.scale,rsep_flat,'o')
#    plt.plot(r/a.scale,pred,'r')
#    plt.ylim(0,2)
#    print(popt)
#







    
        
        
        
        
        