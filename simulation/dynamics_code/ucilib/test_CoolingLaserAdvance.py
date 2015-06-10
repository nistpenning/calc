import CoolingLaserAdvance
import BendKickUpdater
import numpy as np
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array


testCtx = cl.create_some_context(interactive = True)
testQueue = cl.CommandQueue(testCtx)

cla0 = CoolingLaserAdvance.CoolingLaserAdvance(testCtx, testQueue)
cla0.sigma = None
cla1 = CoolingLaserAdvance.CoolingLaserAdvance(testCtx, testQueue)
cla1.sigma = 10.0e-3
cla1.k0 = np.array([0, 0, -2.0 * np.pi / 313.0e-9], dtype = np.float32)

def test_Constructor():
    assert(True)

def test_SetUpdaterParameters():
    global updater
    updater = BendKickUpdater.BendKickUpdater(testCtx, testQueue)
    updater.trapConfiguration.Bz = 1.0e-6


accelerations = []
def test_setAccelerations():
    global accelerations
    accelerations = [cla0, cla1]


def test_initPtcls():
    global ptcls
    ptcls = Ptcls.Ptcls()
    ptcls.numPtcls = 10
    ptcls.init_ptcls()
    # initial mean velocity of particles
    ptcls.ptclList[5] = 20.0 * np.ones_like(ptcls.ptclList[5])

def take_steps(n, dt, ptcls, updater, accelerations):
    xd = cl_array.to_device(updater.queue, ptcls.x(), async = True)
    yd = cl_array.to_device(updater.queue, ptcls.y(), async = True)
    zd = cl_array.to_device(updater.queue, ptcls.z(), async = True)
    vxd = cl_array.to_device(updater.queue, ptcls.vx(), async = True)
    vyd = cl_array.to_device(updater.queue, ptcls.vy(), async = True)
    vzd = cl_array.to_device(updater.queue, ptcls.vz(), async = True)
    qd = cl_array.to_device(updater.queue, ptcls.q(), async = True)
    md = cl_array.to_device(updater.queue, ptcls.m())


    t = 0.0
    for i in range(n):
        updater.update(xd, yd, zd, vxd, vyd, vzd, qd, md,
                accelerations, t, dt, 1)
        t += dt

    xd.get(updater.queue, ptcls.x(), async = True)
    yd.get(updater.queue, ptcls.y(), async = True)
    zd.get(updater.queue, ptcls.z(), async = True)
    vxd.get(updater.queue, ptcls.vx(), async = True)
    vyd.get(updater.queue, ptcls.vy(), async = True)
    vzd.get(updater.queue, ptcls.vz())


def test_takeSteps():
    numSteps = 500
    history = np.ndarray([numSteps, ptcls.numPtcls])
    for i in range(numSteps):
        history[i] = ptcls.vz().copy()
        take_steps(10, 1.0e-6, ptcls, updater, accelerations)
    np.savetxt('vzOfT.dat', history)


# vi: ts=4 sw=4
