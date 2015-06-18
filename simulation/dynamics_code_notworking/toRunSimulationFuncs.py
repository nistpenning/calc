# Function for running simulations

import sys
# sys.path.append('C:\Users\ACKWinDesk\Documents\GitHub\ultracold-ions')
import ucilib.Sim as Sim
import ucilib.TrapConfiguration as TrapConfig
import numpy as np
import ucilib.CoolingLaserAdvance as coolingLsr
import matplotlib.pyplot as plt
import pyopencl as cl
import pyopencl.array as cl_array

import mode_analysis_code as mode
# from time import time
import time


def loadCrystal(filename):
    t = TrapConfig.TrapConfiguration()
    t.Bz = 4.4584
    # V0 in V/m^2
    # t.kz = 2.0 * 1.167e6
    t.kz = 2 * 1.165286288736658e6  # use my value just in case
    delta = 0.0036
    t.kx = -(0.5 + delta) * t.kz
    t.ky = -(0.5 - delta) * t.kz
    # t.theta = 0
    t.theta = 6.266067630465074

    # t.omega = 2.0 * np.pi * 64.0e3
    # t.omega = 2.0 * np.pi * 44.0e3
    t.omega = 2.0 * np.pi * 54.5e3
    # t.omega = 2.0 * np.pi * 50.0e3  #Test High Frequency
    fundcharge = 1.602176565e-19
    # ionMass = 8.9465 * 1.673e-27
    ionMass = 9.012182 * 1.660538921e-27;
    s = Sim.Sim()
    # s.ptcls.set_nptcls(127)
    s.ptcls.set_nptcls(19)
    # s.ptcls.set_nptcls(7)
    s.ptcls.rmax = 1.0e-4
    # s.ptcls.init_ptcls(charge = fundcharge, mass = ionMass) # why is this here if we are loading a file anyway?
    axialFrictionCoeff = None
    angularFrictionCoeff = None
    s.ptcls.ptclList = np.loadtxt(filename, dtype=np.float32)
    s.init_sim(t, axialFrictionCoeff, angularFrictionCoeff,
               recoilVelocity=0.01, scatterRate=1.0e6)

    ##    if s.ptcls.ptclList.shape[0] < 8:
    ##        s.ptcls.ptclList = np.append(s.ptcls.ptclList, fundcharge*np.ones((1,s.ptcls.ptclList.shape[1])),axis=0)
    ##        s.ptcls.ptclList = np.append(s.ptcls.ptclList, ionMass*np.ones((1,s.ptcls.ptclList.shape[1])),axis=0)
    ##    print s.ptcls.ptclList.shape
    s.t = 0.0
    return s


# Need to add cooling laser parameters
def saveParams(sim, dt, nSteps, innerSteps, InitialCrystal, scale, filename):
    t = sim.trapConfiguration
    delta = (t.ky - t.kx) / (2.0 * t.kz)  # reconstruct delta
    N = sim.ptcls.numPtcls  # number of particles
    w = t.omega / (2.0 * np.pi * 1000)  # in kHz
    theta_0 = sim.trapConfiguration.theta
    params = np.array([N, w, delta, theta_0, nSteps, innerSteps, dt, scale])
    np.savetxt(filename, params)
    # with open(filename, "a") as myfile:
    #    myfile.write(InitialCrystal)


# For long simulations, data is too big. Only save positions
def savePtcls_light(sim, filename):
    np.savetxt(filename, sim.ptcls.ptclList[(0, 1, 2, 3, 4, 5), :])  # get everything except mass and charge
    # np.savetxt(filename, sim.ptcls.ptclList[(0,1,2,5),:]) # get positions and velocities
    # np.savetxt(filename, sim.ptcls.ptclList[:3,:(sim.ptcls.numPtcls)])
    # np.savetxt(filename, sim.ptcls.ptclList[:2,:(sim.ptcls.numPtcls)])  # just save plane positions


def saveAfterSim(simdata, filename):
    with file(filename, 'w') as outfile:
        for slice_2d in simdata:
            np.savetxt(outfile, slice_2d)


def savePtcls(sim, filename):
    np.savetxt(filename, sim.ptcls.ptclList)


def buildFileName(path, basename, suffix):
    return path + basename + suffix


def displayCrystal(sim):
    plt.scatter(sim.ptcls.x(), sim.ptcls.y())
    plt.show()


# plt.draw()

def rotate(xy, phase):
    return np.array([
        xy[0] * np.cos(phase) - xy[1] * np.sin(phase),
        xy[0] * np.sin(phase) + xy[1] * np.cos(phase)
    ])


def scrambleZCoords(sim, scale):
    nPtcls = sim.ptcls.ptclList[0].size
    sim.ptcls.ptclList[2] = np.random.standard_normal(nPtcls) * scale


def scrambleXCoords(sim, scale):
    nPtcls = sim.ptcls.ptclList[0].size
    sim.ptcls.ptclList[0] = sim.ptcls.ptclList[0] + np.random.standard_normal(nPtcls) * scale


def scrambleYCoords(sim, scale):
    nPtcls = sim.ptcls.ptclList[0].size
    sim.ptcls.ptclList[1] = sim.ptcls.ptclList[1] + np.random.standard_normal(nPtcls) * scale


def runSim(s, nSteps, innerSteps, dt):
    # axd = cl_array.zeros(s.queue, s.ptcls.x().shape, s.ptcls.x().dtype)
    # ayd = cl_array.zeros(s.queue, s.ptcls.x().shape, s.ptcls.x().dtype)
    # azd = cl_array.zeros(s.queue, s.ptcls.x().shape, s.ptcls.x().dtype)
    for i in range(nSteps):
        # s.take_steps_2(dt, innerSteps,axd, ayd, azd)
        s.take_steps(dt, innerSteps)

    DataLocation = 'D:/PenningSimulationData/'
    date = '2014_5_5_SmallCrystalModes/'
    np.savetxt(DataLocation + date + 'final_copy.dat', s.ptcls.ptclList)


# #        if not (i % 10):
# #            ti = time.time()
# #            print "i: ", i, "wall time: ", ti - t0, "sim. time: ", s.t
# #            log[i/1000,0] = ti - t0
# #            log[i/1000,1] = s.t

def MakeModeCodeCrystal(N=127, Vtrap=(0.0, -1750.0, -2000.0), Ctrap=1.0,
                        ionmass=None, B=4.4588, frot=180., Vwall=1., wall_order=2,
                        quiet=True):
    """
    #TODO work on obtaining the vtrap, B, etc from a trapconfiguration class.
    All of the parameters below are the parameters the mode code wants for a crystal.

    :param N:
    :param Vtrap:
    :param Ctrap:
    :param ionmass:
    :param B:
    :param frot:
    :param Vwall:
    :param wall_order:
    :param quiet:
    :return:
    """
    a = mode.ModeAnalysis(N=N, Vtrap=Vtrap, Ctrap=Ctrap, ionmass=ionmass, B=B, frot=frot,
                          Vwall=Vwall, wall_order=wall_order)
    a.run()
    return (a)


def GetXYFromModeCrystal(crystal):
    N=crystal.Nion
    x =crystal.u[:N]
    x=np.array(x,dtype=np.float32)
    x = x.reshape(1, crystal.Nion)

    y =crystal.u[N:]
    y=np.array(y,dtype=np.float32)
    y = y.reshape(1, crystal.Nion)

    z = np.zeros((1, crystal.Nion), dtype= np.float32)
    v = np.zeros((3, crystal.Nion), dtype = np.float32)
    q = crystal.q * (np.ones((1, crystal.Nion), dtype=np.float32))
    m = crystal.m_Be * (np.ones((1, crystal.Nion), dtype= np.float32))


    particles = np.concatenate((x, y, z, v, m, q), axis=0)
    return particles


def GetEvalsFromModeCrystal(crystal, dimension="axial", experimentalunits=False):
    dimension = dimension.lower()
    if experimentalunits is False and dimension[0] is "a":
        return crystal.axialEvals
    if experimentalunits is False and dimension[0] is "p":
        return crystal.planarEvals
    if experimentalunits is True and dimension[0] is "a":
        return crystal.axialEvalsE
    if experimentalunits is True and dimension[0] is "p":
        return crystal.planarEvalsE


def GetEvectsFromModeCrystal(crystal, dimension="axial"):
    dimension = dimension.lower()
    if dimension[0] is "a":
        return crystal.axialEvects
    if dimension[0] is "planar":
        return crystal.planarEvects
