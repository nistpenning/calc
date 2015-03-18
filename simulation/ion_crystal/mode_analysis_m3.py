# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 14:12:58 2014

@author: jgb
"""
from __future__ import division
from scipy.constants import pi
import numpy as np
from scicons import m_Be, q, k_e
import scipy.optimize as optimize
import matplotlib.pyplot as plt

class ModeAnalysis():
    def __init__(self, shells=4, Vtrap = [0.0,-200.0,-400.0], Ctrap = 1.0, 
                 fz=1000, B=4.4609, frot=60, Vwall=1, wall_order=2, mult=1e14):
        self.quiet = False        
        self.Evects = 0 #Eigen vectors
        self.Evals = 0 #Eigen frequencies
        self.dens = 0
        self.avg = 0
        
        self.shells = shells
        self.Nion = 1 + 6* np.sum(range(1,shells+1))
        self.u0 = np.empty(2*self.Nion) #for array of ion positions first half is x, last is y
        self.u = np.empty(2*self.Nion) #for array of ion positions first half is x, last is y
        self.scale = 0
        
        #trap definitions
        self.B = B
        self.m = m_Be * np.ones(self.Nion)
        self.wcyc = q*B/m_Be
        self.C = Ctrap * np.array([[0.0756,0.5157,0.4087],
                   [-0.0001,-0.005,0.005],
                   [1.9197e3,3.7467e3,-5.6663e3],
                   [0.6738e7,-5.3148e7,4.641e7]]) #axial trap coefficients
        self.relec = 0.01 #rot wall electrode distance in meters
        self.Vtrap = np.array(Vtrap) #Vend, Vmid, Vcenter
        self.Coeff = np.dot(self.C,self.Vtrap)
        
        self.wz = np.sqrt(2*q*self.Coeff[2]/m_Be)
        self.wrot = 2*pi*frot*1e3
        self.wmag = 0.5*(self.wcyc-np.sqrt(self.wcyc**2-2*self.wz**2))

        self.V0 = (0.5*m_Be*self.wz**2)/q
        self.Vw = self.V0 *0.045*Vwall/1000 #V_T divides in MATLAB
        self.mult = mult
        
        #wall order
        if wall_order == 2:
            self.Cw2 = q*Vwall*1612
            self.Cw3 = 0
        if wall_order == 3:
            self.Cw2 = 0
            self.Cw3 = q*Vwall*3e4

    
    def run(self):
        mins = 10e-6
        res = 0.1e-6
        if self.wmag > self.wrot:
            print("Warning: Rotation frequency below magnetron frequency of {0.1f}".format(self.wmag))
        self.u0[:] = self.find_scaled_lattice_guess(mins, res)
        #self.show_crystal(self.u0)
        print "Calculating Equilibrium Positions"
        self.u = self.find_eq_pos(self.u0)
        #self.show_crystal(self.u)
        print "Calculate transverse axial modes"
        self.Evals,self.Evects = self.calc_axial_modes(self.u)
        
        #sort arrays
        sort_ind = np.argsort(np.abs(self.Evals))[::-1]
        self.Evals = self.Evals[sort_ind]/(2*pi*1e3) #units of kHz
        self.Evects = self.Evects[:,sort_ind]
        
        #plt.plot((self.Evals,self.Evals),(np.ones(np.size(self.Evals)),np.zeros(np.size(self.Evals))),linestyle="-",color='black')
        #plt.axis([800,1010,0,1.1])
        
        self.show_crystal_modes(self.u, self.Evects,3)
        
    def run_quiet(self):
        self.quiet = True
        mins = 10e-6
        res = 0.1e-6
        if self.wmag > self.wrot:
            print("Warning: Rotation frequency below magnetron frequency of {0:.1f}".format(float(self.wmag/2*pi)))
            return 0
        self.u0[:] = self.find_scaled_lattice_guess(mins, res)
        self.u = self.find_eq_pos(self.u0)
        
        self.Evals,self.Evects = self.calc_axial_modes(self.u)
        
        #sort arrays
        sort_ind = np.argsort(np.abs(self.Evals))[::-1]
        self.Evals = self.Evals[sort_ind]/(2*pi*1e3) #units of kHz
        self.Evects = self.Evects[:,sort_ind]
        self.r,dx,dy,self.rsep = self.find_radial_separation(self.u)
       

    def find_scaled_lattice_guess(self, mins, res):
        uthen = self.generate_2D_hex_lattice(self.shells, mins)
        pthen = self.pot_energy(uthen)
        
        for scale in np.arange(mins+res, (mins/res)*mins, res):
            uguess = self.generate_2D_hex_lattice(self.shells, scale)
            pnow = self.pot_energy(uguess)
            if pnow >= pthen:
                #print "find_scaled_lattice: Minimum found"
                #print "initial scale guess: " + str(scale)
                self.scale = scale
                return uthen
            uthen = uguess
            pthen = pnow
        self.scale = scale
        #print "find_scaled_lattice: no minimum found, returning last guess"
        return uthen
            
    def generate_2D_hex_lattice(self, shells = 1, scale = 1):
        posvect =  np.array([0.0,0.0])#always a point [0,0]
        
        for s in range(1,self.shells+1):
            posvect = np.append(posvect, self.add_hex_shell(s))
        posvect = scale * posvect
    
        posvect = np.transpose(np.reshape(posvect, (self.Nion,2)))
        return np.reshape(posvect,2*self.Nion)
    
    def add_hex_shell(self, s):
        a = range(s,-s-1,-1)
        a.extend(-s*np.ones(s-1))
        a.extend(range(-s,s+1))
        a.extend(s*np.ones(s-1))
    
        b = range(0,s+1)
        b.extend(s*np.ones(s-1))
        b.extend(range(s,-s-1,-1))
        b.extend(-s*np.ones(s-1))
        b.extend(range(-s,0))
        
        x = np.sqrt(3)/2.0 * np.array(b)
        y = 0.5 * np.array(b) + np.array(a)
        
        pair = np.column_stack((x,y)).flatten()
        return pair
    
    def pot_energy(self, pos_array):
        w = self.wrot
        m = self.m[0]
        N = int(pos_array.size/2)
        
        x = pos_array[0:N]
        y = pos_array[N:]
        
        dx = x.reshape((x.size,1)) - x
        dy = y.reshape((y.size,1)) - y
        rsep = np.sqrt(dx**2 + dy**2)   
        
        #dxsq = np.array([(i-x)**2 for i in x])
        #dysq = np.array([(i-y)**2 for i in y])

        #dxsq = [(i-x)**2 for i in x]
        #dysq = [(i-y)**2 for i in y]
        
        #rsep = np.sqrt(dxsq + dysq)        
        Vc = 1/rsep
        Vc[np.isinf(Vc)|np.isnan(Vc)] = 0
        
        V = 0.5*(-m*w**2 - q*self.Coeff[2] + q*self.B*w) * np.sum((x**2 + y**2)) \
                   - q*self.Coeff[3] *  np.sum((x**2 + y**2)**2) \
                   + np.sum( self.Cw2 *(x**2 - y**2)) \
                   + np.sum( self.Cw3 *(x**3 - 3*x*y**2)) \
                   + 0.5 * k_e * q**2 *np.sum(Vc)
        
        return self.mult*V
        
    
    def force_penning(self, pos_array):
        w = self.wrot
        m = self.m[0]
        N = int(pos_array.size/2)
        
        x = pos_array[0:N]
        y = pos_array[N:]
        
        dx = x.reshape((x.size,1)) - x
        dy = y.reshape((y.size,1)) - y
        rsep = np.sqrt(dx**2 + dy**2)   
                            
        #Calculate coulomb force on each ion
        Fc = rsep**(-2)
        Fc[np.isinf(Fc)|np.isnan(Fc)] = 0
        fx = np.float64((dx/rsep)*Fc) 
        fx[np.isinf(fx)|np.isnan(fx)] = 0
        fy  = np.float64((dy/rsep)*Fc)
        fy[np.isinf(fy)|np.isnan(fy)] = 0
        
        #total force on each ion
#       Ftrapx = (m*w**2 + q*self.V0 - 2*q*self.Vw - q*self.B* w) * x
#       Ftrapy = (m*w**2 + q*self.V0 + 2*q*self.Vw - q*self.B* w) * y
        Ftrapx = (-m*w**2 -q*self.Coeff[2] +q*self.B*w +2*self.Cw2) * x \
                    - 4*q*self.Coeff[3]*(x**3+x*y**2) + 3*self.Cw3*(x**2-y**2)
        Ftrapy = (-m*w**2 -q*self.Coeff[2] +q*self.B*w -2*self.Cw2) * y \
                    - 4*q*self.Coeff[3]*(y**3+y*x**2) - 6*self.Cw3*x*y

        #Ftrap =  (m*w**2 + q*self.V0 - 2*q*self.Vw - q*self.B* w) * pos_array
        Fx = -(k_e*q**2)*np.sum(fx, axis = 1) + Ftrapx
        Fy = -(k_e*q**2)*np.sum(fy, axis = 1) + Ftrapy
        
        Fx[np.abs(Fx) < 1e-24] = 0
        Fy[np.abs(Fy) < 1e-24] = 0

        Fx = self.mult*Fx
        Fy = self.mult*Fy
        #self.F = np.array([Fx,Fy]).flatten()
        
        return np.array([Fx,Fy]).flatten()

    def hessian_penning(self, pos_array):
        w = self.wrot
        m = self.m[0]        
        N = int(pos_array.size/2)
        H = np.empty((2*N,2*N))
        
        x = pos_array[0:N]
        y = pos_array[N:]
        
        dx = x.reshape((x.size,1)) - x
        dy = y.reshape((y.size,1)) - y
        rsep = np.sqrt(dx**2 + dy**2)
        rsep5 = self.nan_to_zero(rsep**-5)
        
        dxsq = dx**2
        dysq = dy**2
        
        #in these calcs, diagnonal indices should always be zero
        #np.fill_diagonal(dxsq, 0)
        #np.fill_diagnoal(dysq, 0)
        
        H1 = np.mat((k_e*q**2) * (rsep**2 - 3*dxsq) * rsep5)
        
        H5 = np.mat(np.diag((-m*w**2 -q*self.Coeff[2] +q*self.B*w +2*self.Cw2) \
                        + 6*self.Cw3*x - 4*q*self.Coeff[3]*(3*x**2+y**2) \
                        - k_e*q**2 * np.sum((rsep**2 - 3*dxsq) * rsep5, axis=0)))        
                        
        H2 = np.mat(-(k_e*q**2) * (rsep**2 - 3*dysq) * rsep5)
        H6 = np.mat(np.diag((-m*w**2 -q*self.Coeff[2] +q*self.B*w -2*self.Cw2) \
                        - 6*self.Cw3*x - 4*q*self.Coeff[3]*(3*y**2+x**2) \
                        - k_e*q**2 * np.sum((rsep**2 - 3*dysq) * rsep5, axis=0)))
        
        H3 = np.mat(-3*k_e*q**2 * dx * dy * rsep5)
        H4 = np.mat(np.diag(3*k_e*q**2 * np.sum(dx * dy * rsep5, axis=0)))
        
        H11 = H1 + H5
        H22 = H2 + H6
        H12 = H3 + H4
                        
        H = np.bmat([[H11, H12],[H12, H22]])
        
        H = np.asarray(H)
        
        #H[np.abs(H) < 1e-20] = 0
                
        return self.mult*H
        
    def find_eq_pos(self, u0):
        fun_tolerance = 1e-37
        
        #return optimize.fmin_bfgs(self.pot_energy, u0, maxiter = 1000)
        
        #out = optimize.fmin_ncg(self.pot_energy, u0, self.force_penning, fhess=self.hessian_penning, 
        #                        maxiter=1000, ftol=fun_tolerance)
        
        #guessadd = np.array([4e-6*np.random.random() for i in u0])
                                
#        out = optimize.minimize(self.pot_energy, u0, method='Newton-CG', jac=self.force_penning, 
#                                hess=self.hessian_penning, tol=fun_tolerance, options={'xtol': 1e-37,'disp': True})

        out = optimize.minimize(self.pot_energy, u0, method='BFGS', jac=self.force_penning,
                                options={'gtol': fun_tolerance, 'disp': not self.quiet})
        
        return out.x
                
    def calc_axial_modes(self, pos_array):
        m = self.m[0]
        N = int(pos_array.size/2)
        A = np.empty((N,N))
        
        x = pos_array[0:N]
        y = pos_array[N:]
        
        dx = x.reshape((x.size,1)) - x
        dy = y.reshape((y.size,1)) - y
        rsep = np.sqrt(dx**2 + dy**2)
        rsep3 = self.nan_to_zero(rsep**-3)
        
        A1 = np.diag((2*q*self.Coeff[2] - k_e*q**2*np.sum(rsep3, axis=0)))
        A2 = k_e*q**2 * rsep3
        A[:] = (A1 + A2)/m
        
        Eval,Evect = np.linalg.eig(A)
        
        Eval = np.lib.scimath.sqrt(Eval)
        
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
        my_array[np.isinf(my_array)|np.isnan(my_array)] = 0
        return my_array
        
    def show_crystal(self, pos_vect):
        plt.plot(1e6*pos_vect[0:self.Nion],1e6*pos_vect[self.Nion:],'.')
        #ax.set_aspect('equal')
        plt.xlabel('x position [um]')
        plt.ylabel('y position [um]')
        plt.axes().set_aspect('equal')
        plt.axis([-200,200,-200,200])     
        
        plt.show()
        
    def show_crystal_modes(self, pos_vect, Evects, modes):
        plt.figure(1)
        
        for i in range(modes):
            plt.subplot(1,modes,i+1,aspect='equal')
            plt.scatter(1e6*pos_vect[0:self.Nion],1e6*pos_vect[self.Nion:],c=Evects[:,i], vmin=-.25, vmax=0.25)
            #ax.set_aspect('equal')
            plt.xlabel('x position [um]')
            plt.ylabel('y position [um]')
            plt.axis([-200,200,-200,200])
            #plt.axes().set_aspect('equal')
            
        plt.tight_layout()
        
        
    
    def save_positions(self, u):
        np.savetxt("py_u.csv", u, delimiter = ",")
        
    def crystal_spacing_fit(self, r, offset, curvature):
        return np.sqrt(2/(np.sqrt(3)*offset*np.sqrt(1-(r*curvature)**(2))))
                
if __name__ == "__main__":
    a = ModeAnalysis(shells=7 ,Vtrap=[-0.0,-495.0,-596.0], Ctrap = 1.0, frot=58.0, Vwall= 0.10, wall_order=2)
    a.run_quiet()
    Evals = a.Evals
    Evect = a.Evects
#    plt.close()
#    plt.plot((Evals,Evals),(np.ones(np.size(Evals)),np.zeros(np.size(Evals))),linestyle="-",color='black')
#    plt.axis([800,900,0,1.1])
#    plt.ylabel('arb')
#    plt.xlabel('Mode Frequency [kHz]')
#    a.show_crystal_modes(a.u, a.Evects, 3)
#    print(np.max(Evals))
#    print(np.min(Evals))
#    plt.close()
    a.show_crystal(a.u)
#    print(a.Nion)
    
#   get the data    
    r = np.zeros(np.shape(a.r))
    rsep = np.zeros(np.shape(a.rsep))
    r[:] = a.r
    rsep[:]= np.transpose(a.rsep)/a.scale

#   process into a single 1D array -- use braodcasting to make a r into repeating 
#   2D array (faster than a loop)
    r_full = r + np.zeros(np.shape(rsep)) 
    rsep_flat = rsep.flatten()
    r_flat = r_full.flatten()
    
#   create mask for points
    rsep_flat[rsep_flat>1.4] = np.nan
    rsep_flat[rsep_flat<0.1] = np.nan
    mask = np.isnan(rsep_flat)
    r_fit = r_flat[~mask]
    rsep_fit = rsep_flat[~mask]
    
    popt,pcov = optimize.curve_fit(a.crystal_spacing_fit, r_fit, rsep_fit,
                                    p0=[1.25,2000])
    
    pred = a.crystal_spacing_fit(r, *popt)
    
    plt.plot(r_flat/a.scale,rsep_flat,'o')
    plt.plot(r/a.scale,pred,'r')
    plt.ylim(0,2)
    print(popt)








    
        
        
        
        
        