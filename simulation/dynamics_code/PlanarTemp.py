##
## Input one my T=0 crystal with small z positions in normal mode expansion
## want to see evolution of normal modes
## No cooling beams, just trap potentials and B field
## 

import sys
sys.path.append('C:\Users\ACKWinDesk\Documents\GitHub\ultracold-ions')
import os
os.environ['PYOPENCL_CTX'] = '1' # CPU is doing computation
import ucilib.Sim as Sim
import ucilib.TrapConfiguration as TrapConfig
import numpy as np
import ucilib.CoolingLaserAdvance as coolingLsr
import matplotlib.pyplot as plt
from toRunSimulationFuncs import*   # import functions that all simulation runs share
#from time import time
import time

import cProfile


DataLocation = 'D:/PenningSimulationData/'
#date = '2014_7_17_PlanarTemp/'
#date = '2014_9_10_PlanarTemp/'
#date = '2014_10_16_PlanarTemp/'
date = '2014_11_07_PlanarTemp/'


mode = 1
G =  4.5e-5
Vwall = -400
   
#InitialCrystal = '7_64_' + str(Vwall) + '_ZeroTemp.dat'
InitialCrystal = '19_54.5_' + str(Vwall) + '_ZeroTemp.dat'
#InitialCrystal = '19_54.5_-400_Continuation10_06.dat'
#InitialCrystal = '127_44_' + str(Vwall) + '_ZeroTemp.dat'
s = loadCrystal(DataLocation + date + InitialCrystal) # this has theta = 0\
s.updateDelta(-G*Vwall)

scale = 0

s.spin_up()
s.accList = s.accList[0:2]
s.updater.axialFrictionCoeff=None
s.updater.angularFrictionCoeff=None

claAlongZ = coolingLsr.CoolingLaserAdvance(s.ctx, s.queue)
claAlongZ.sigma = None
claAlongZ.k0 = np.array([0, 0, -2.0 * np.pi / 313.0e-9], dtype = np.float32)
s.accList.append(claAlongZ)

perpCooling = coolingLsr.CoolingLaserAdvance(s.ctx, s.queue)
#perpCooling.sigma = 30.0e-6
perpCooling.sigma = 15.0e-6
perpCooling.k0 = np.array([2.0 * np.pi / 313.0e-9, 0, 0], dtype = np.float32)
#perpCooling.k0 = np.array([-2.0 * np.pi / 313.0e-9, 0, 0], dtype = np.float32)
#perpCooling.x0 = np.array([0, 18.0e-6, 0], dtype = np.float32)
perpCooling.x0 = np.array([0, -9.0e-6, 0], dtype = np.float32)
#perpCooling.S = 1
perpCooling.S = 0.3
s.accList.append(perpCooling)

nSteps = 1000000     # Number of snapshots
innerSteps = 100   # Inner number of steps 
dt = 5e-10
 
saveParams(s,dt,nSteps,innerSteps,InitialCrystal,scale,DataLocation+date+"params.dat")  # save run parameters
thetas = np.array([s.trapConfiguration.theta]) 

simData = np.ndarray([nSteps,6,s.ptcls.numPtcls])
log = np.ndarray([nSteps/1,2])

t0 = time.time()
LASERON = 1
binsize = 40000
for i in range(nSteps):
        
        if not (i % 10000):
                ti = time.time()
                print "i: ", i, "wall time: ", ti - t0, "sim. time: ", s.t
                log[i/1000,0] = ti - t0
                log[i/1000,1] = s.t

                        
        simData[i] = s.ptcls.ptclList[(0,1,2,3,4,5),:]  # add current step to history


##        if not (i % binsize):
##                if LASERON == -1: # if laser is off, turn it on
##                        s.accList.append(claAlongZ)
##                        s.accList.append(perpCooling)
##                if LASERON == 1: # if laser is on
##                        s.accList = s.accList[0:2]
##                LASERON = -1*LASERON;

        if i == (nSteps-1):
                cProfile.run('s.take_steps(dt, innerSteps)','profilestats')
        else:
                s.take_steps(dt, innerSteps)
        thetas = np.append(thetas, s.trapConfiguration.theta)


np.savetxt(DataLocation+date+'thetas.dat', thetas) #save thetas
np.savetxt(DataLocation+date+'log3' + str(mode) + '.log', log) #save thetas

saveAfterSim(simData, DataLocation+date+ 'Data.dat')
#saveAfterSim(simData, DataLocation+date+ 'Data' + str(mode) + '.dat')
#saveAfterSim(simData, DataLocation+date+ 'Data' + str(ind+1) + '.dat')
print "elapsed wall time: ", time.time()-t0
