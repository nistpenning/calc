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
date = '2014_5_13_SmallCrystalModes_2/'
#date = '2014_5_12_LargeCrystalModes/'
#date = '2014_5_4_ProfilingTest/'

for mode in range(1,8):

    #InitialCrystal = '19_54.5_-80_Amode' + str(mode) + '.dat'
    #InitialCrystal = '127_44_-80_Amode' + str(mode) + '.dat'
    InitialCrystal = '7_64_-160_Amode' + str(mode) + '.dat'
    s = loadCrystal(DataLocation + date + InitialCrystal) # this has theta = 0\

    #InitialCrystal = '19_54.5_-80_eq.dat'
    #InitialCrystal = '127_44_-80_truemin.dat'
    #s = loadCrystal('C:/Users/ACKWinDesk/Desktop/PenningSimulationData/MatlabData/' + InitialCrystal) # this has theta = 0\
    

    scale = 5e-6

    s.spin_up()
    s.accList = s.accList[0:2]
    s.updater.axialFrictionCoeff=None
    s.updater.angularFrictionCoeff=None

    nSteps = 100000     # Number of snapshots
    innerSteps = 100   # Inner number of steps 
    dt = 5e-9
     
    saveParams(s,dt,nSteps,innerSteps,InitialCrystal,scale,DataLocation+date+"params.dat")  # save run parameters
    thetas = np.array([s.trapConfiguration.theta]) 

    simData = np.ndarray([nSteps,6,s.ptcls.numPtcls])
    log = np.ndarray([nSteps/1,2])
    #sys.exit()
    t0 = time.time()

##    cProfile.run('runSim(s,nSteps,innerSteps,dt)','restats')
##    import pstats
##    p = pstats.Stats('restats')
##    p.sort_stats('cumulative').print_stats(30)
    
    for i in range(nSteps):
        
        if not (i % 1000):
            ti = time.time()
            print "i: ", i, "wall time: ", ti - t0, "sim. time: ", s.t
            log[i/1000,0] = ti - t0
            log[i/1000,1] = s.t

        simData[i] = s.ptcls.ptclList[(0,1,2,3,4,5),:]  # add current step to history
        #cProfile.run('s.take_steps(dt, innerSteps)')
        s.take_steps(dt, innerSteps)
        thetas = np.append(thetas, s.trapConfiguration.theta)


    np.savetxt(DataLocation+date+'thetas.dat', thetas) #save thetas
    np.savetxt(DataLocation+date+'log' + str(mode) + '.log', log) #save thetas

    saveAfterSim(simData, DataLocation+date+ 'Data' + str(mode) + '.dat')
    print "elapsed wall time: ", time.time()-t0
