# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 19:36:34 2015

@author: ACKWinDesk

These are methods that should be in a total 
analysis package for the simulation code. Extracts mode energies
and does a fourier analysis. Still untested as having trouble installing
new simulation files
"""

import mode_analysis_code as MC
import numpy as np

def readFromFile():
    """Convert file to pctl list"""
    pass

def loadData(self):
    """Read all data from filename"""
    self.readFromFile()
    pass    

def fourierAnalysis():
    """Load all data and compute power spectral densities"""
    self.loadData()
    
def kineticEnergy():
    """Calculate kinetic energy of all particles"""
    pass

def potentialEnergy():
    """Calculate potential energy of all particles"""
    pass
    
def exciteMode():
    pass


def spin_down(u, v, w):
    """Find velocities in rotating frame (move to that frame)"""
    N = int(u.size/2)
    radii = np.sqrt(u[0:N]**2 + u[N:]**2)
    velocities = w*radii
    for i in range(N):
        rot = [-u[i+N], u[i]]
        rot = rot/numpy.linalg.norm(rot)
        # counter velocity to move to rotating frame
        rot = -velocities[i]*rot;
        v[i] = v[i] + rot[0]
        v[i+N] = v[i+N] + rot[1]
    return v  
        
def axialModeBasis(z, Evect):
    """Project axial motion into normal coordinates of axial modes"""
    N = int(u.size/2)
    norm_coords = zeros(1,N);

for m = modes
    norm_coords(m) = dot(z',E(:,m));

def rotate(u, theta):
    """Rotates coordinates by theta"""
    N = int(u.size/2)
    x = u[0:N]
    y = u[N:]

    xnew = x*np.cos(theta) - y*np.sin(theta)
    ynew = x*np.sin(theta) + y*np.cos(theta)
    
    return np.hstack((x0, y0))



aTrap = ModeAnalysis(N=nions, Vtrap=[0.0, -1750.0, -2000.0], Ctrap=1.0,
                     frot=180, Vwall=2, wall_order=2)
aTrap.run()
V_0 = aTrap.pot_energy(aTrap.u)  # initial potential energy of crystal
# need to rotate to get same major axis as dominics code (his is along x-axis)
# before or after solving?
uE = rotate(self.uE, pi/2)

# Run Dominic's code
# get output in some way
#  Let's call us a time series of crystals (each row is a snapshot)
#  Let's call vs a time series of crystal planar velocities
#  Let's call zs a time series of axial ion motion


norm_coords = np.zeros((steps,aTrap.Nion)) 
norm_coords_planar = np.zeros((steps, 4*aTrap.Nion)) 
v = np.zeros((steps, 2*aTrap.Nion))
PlanarMotion = np.zeros((steps, 2*aTrap.Nion))
EnergyPlanarMode = np.zeros((steps, 2*aTrap.Nion))
EnergyAxialMode = np.zeros((steps, aTrap.Nion))


PPE = []  # crystal potential energy at every snapshot 
# for every snapshot in simulation
for s in range(steps):
    vrot[s, :] = spin_down(us[s,:], vs[s,:], w)
    us[s, :] = rotate(us[s,:],-thetas[i])  # rotate crystal back
    #PPE.append(0.5*m*wz^2*l0^2*( aTrap.pot_energy(us[s,:]) - V_0 )) # this needs to be fixed for dimensionless
    norm_coords[s,:] = modeBasis(zs[s,:], aTrap.axialEvects) 
    PlanarMotion[s, :] = us[s, :] - uE # subtract off equilibrium positions  
    qqdot = np.array([PlanarMotion[s, :], vrot[s,:])
    # need to get these dimensions right so matrix product works
    norm_coords_planar[i,:] = numpy.linalg.inv(aTrap.planarEvects)*qqdot 
        


for j in range(2*aTrapNion):
    EnergyPlanarMode[:, j] = 0.5*m*wz^2/kB*l0**2*(np.absolute(norm_coords_planar[:,2*j])**2+
                                                 np.absolute(norm_coords_planar[:,2*j+1])**2)

for j in range(aTrapNion):
    EnergyAxialMode[:, j] = 0.5*m*wz^2/kB*l0**2*(np.absolute(norm_coords_planar[:,2*j])**2+
                                                 np.absolute(norm_coords_planar[:,2*j+1])**2)


# MATLAB PSD CODE
InSteps = 1e2   # Inner steps in simulation
OutSteps = 1e6  # Outer number of steps in simulation
dt = 5e-10      # time step

freq = np.arange(0, 0.5/(InSteps*dt), 1/(InSteps*dt*OutSteps))
freq = (0:0.5/(InSteps*dt)/(params(5)/2):0.5/(InSteps*dt))
freq = 1.0e-6*freq(1:end-1); % chop of last one, because of Matlab ranges...

# Calculate PSD for Axial Motion
spectra = abs(fft(zs)).^2;
Zpsd = sum(spectra, 2);
Zpsd = Zpsd(1:(length(Zpsd)/2))+Zpsd(end:-1:length(Zpsd)/2+1);

# Calculate PSD for Planar Motion
motion = us - repmat(us(1,:),params(5),1); % subtract off equilibrium positions
spectra = abs(fft(motion)).^2;
Ppsd = sum(spectra, 2);
Ppsd = Ppsd(1:(length(Ppsd)/2))+Ppsd(end:-1:length(Ppsd)/2+1);

    
