import BendKickUpdater
import TrapAcc
import CoulombAcc
import Ptcls
import TrapConfiguration
import FrictionAcc
import HeatingAcc

import numpy
import pyopencl as cl
import pyopencl.array as cl_array
import copy
import cProfile

class Sim():

    def __init__(self, ctx = None, queue = None):
        self.accList = []
        self.ptcls = Ptcls.Ptcls()
        self.t = 0.0

        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties=cl.command_queue_properties.PROFILING_ENABLE)

        self.trapConfiguration = TrapConfiguration.TrapConfiguration()
        self.updater = BendKickUpdater.BendKickUpdater(self.ctx, self.queue)
        self.updater.trapConfiguration = self.trapConfiguration

    def copy(self):
        theCopy = Sim(self.ctx, self.queue)
        theCopy.accList = self.accList
        theCopy.updater = self.updater
        theCopy.t = copy.deepcopy(self.t)
        theCopy.trapConfiguration = copy.deepcopy(self.trapConfiguration)
        theCopy.updater.trapConfiguration = theCopy.trapConfiguration
        theCopy.ptcls.ptclList = copy.deepcopy(self.ptcls.ptclList)
        theCopy.init_sim(theCopy.trapConfiguration,
                theCopy.updater.axialDampingCoefficient,
                theCopy.updater.angularDampingCoefficient)
        return theCopy
    
    def updateDelta(self,newdelta):
        self.trapConfiguration.changeDelta(newdelta)
        self.updater.trapConfiguration.changeDelta(newdelta)
    
    def init_sim(self, trapConfig, axialDampingCoefficient,
            angularDampingCoefficient, recoilVelocity = None,
            scatterRate = None):

        if recoilVelocity == None:
            recoilVelocity = 0.1
        if scatterRate == None:
            scatterRate = 0.0

        self.trapConfiguration = trapConfig
        self.updater.trapConfiguration = trapConfig
        self.updater.axialDampingCoefficient = axialDampingCoefficient
        self.updater.angularDampingCoefficient = angularDampingCoefficient
        self.recoilVelocity = recoilVelocity
        self.scatterRate = scatterRate

        self.accList = []

        self.accList.append(CoulombAcc.CoulombAcc(self.ctx, self.queue))

        self.trapAcc = TrapAcc.TrapAcc(self.ctx, self.queue)
        self.trapAcc.trapConfiguration = trapConfig
        self.accList.append(self.trapAcc)

        # Make accelerations arrays for BendKickUpdater (load ptcls before calling init_sim)
        self.updater.SetUpAs(self.ptcls.x().shape, self.ptcls.x().dtype)
        
        self.t = 0.0

    def take_steps(self, dt, numSteps = 1):
        xd = cl_array.to_device(self.queue, self.ptcls.x(), async = True)
        yd = cl_array.to_device(self.queue, self.ptcls.y(), async = True)
        zd = cl_array.to_device(self.queue, self.ptcls.z(), async = True)
        vxd = cl_array.to_device(self.queue, self.ptcls.vx(), async = True)
        vyd = cl_array.to_device(self.queue, self.ptcls.vy(), async = True)
        vzd = cl_array.to_device(self.queue, self.ptcls.vz(), async = True)
        qd = cl_array.to_device(self.queue, self.ptcls.q(), async = True)
        md = cl_array.to_device(self.queue, self.ptcls.m())
        #cProfile.run('self.updater.update(xd, yd, zd, vxd, vyd, vzd, qd, md,self.accList, self.t, dt, numSteps)')
        self.t = self.updater.update(xd, yd, zd, vxd, vyd, vzd, qd, md,
                self.accList, self.t, dt, numSteps)

        self.queue.finish()
        xd.get(self.queue, self.ptcls.x(), async = True)
        yd.get(self.queue, self.ptcls.y(), async = True)
        zd.get(self.queue, self.ptcls.z(), async = True)
        vxd.get(self.queue, self.ptcls.vx(), async = True)
        vyd.get(self.queue, self.ptcls.vy(), async = True)
        vzd.get(self.queue, self.ptcls.vz())

    def take_steps_2(self, dt, numSteps ,axd, ayd, azd):
        xd = cl_array.to_device(self.queue, self.ptcls.x(), async = True)
        yd = cl_array.to_device(self.queue, self.ptcls.y(), async = True)
        zd = cl_array.to_device(self.queue, self.ptcls.z(), async = True)
        vxd = cl_array.to_device(self.queue, self.ptcls.vx(), async = True)
        vyd = cl_array.to_device(self.queue, self.ptcls.vy(), async = True)
        vzd = cl_array.to_device(self.queue, self.ptcls.vz(), async = True)
        qd = cl_array.to_device(self.queue, self.ptcls.q(), async = True)
        md = cl_array.to_device(self.queue, self.ptcls.m())
        #cProfile.run('self.updater.update(xd, yd, zd, vxd, vyd, vzd, qd, md,self.accList, self.t, dt, numSteps)')
        self.t = self.updater.update(xd, yd, zd, vxd, vyd, vzd, qd, md,
                self.accList, self.t, dt, numSteps,axd,ayd,azd)

        self.queue.finish()
        xd.get(self.queue, self.ptcls.x(), async = True)
        yd.get(self.queue, self.ptcls.y(), async = True)
        zd.get(self.queue, self.ptcls.z(), async = True)
        vxd.get(self.queue, self.ptcls.vx(), async = True)
        vyd.get(self.queue, self.ptcls.vy(), async = True)
        vzd.get(self.queue, self.ptcls.vz())

    def spin_up(self):
        radii = numpy.sqrt(self.ptcls.x() ** 2 + self.ptcls.y() ** 2)
        velocities = self.updater.trapConfiguration.omega * radii
        for i in range(0, self.ptcls.ptclList.shape[1]):
            v = numpy.array([-self.ptcls.y()[i], self.ptcls.x()[i]])
            v = v / numpy.linalg.norm(v)
            v = velocities[i] * v
            self.ptcls.vx()[i] = v[0] +  self.ptcls.vx()[i] # move velocities to rotating frame
            self.ptcls.vy()[i] = v[1] +  self.ptcls.vy()[i]  

    def radial_velocities(self):
        radVel = []
        for i in range(0, self.ptcls.ptclList.shape[1]):
            radialUnitVec = numpy.array([self.ptcls.x()[i], self.ptcls.y()[i]])
            radialUnitVec = radialUnitVec / numpy.linalg.norm(radialUnitVec)
            radVel.append(numpy.inner(radialUnitVec, numpy.array([self.ptcls.vx()[i], self.ptcls.vy()[i]])))
        return numpy.array(radVel)

    def angular_velocities(self):
        angVel = []
        for i in range(0, self.ptcls.ptclList.shape[1]):
            angularUnitVec = numpy.array([-self.ptcls.y()[i], self.ptcls.x()[i]])
            angularUnitVec = angularUnitVec / numpy.linalg.norm(angularUnitVec)
            angVel.append(numpy.inner(angularUnitVec, numpy.array([self.ptcls.vx()[i], self.ptcls.vy()[i]])))
        return numpy.array(angVel)

    def radii(self):
        return numpy.array([numpy.linalg.norm(numpy.array([self.ptcls.x()[i], self.ptcls.y()[i]])) for i in range(0, self.ptcls.ptclList.shape[1])])

    def xyInRotatingFrame(self):
        theta = self.updater.trapConfiguration.theta
        return numpy.array(
                [numpy.cos(theta) * self.ptcls.x() + numpy.sin(theta) * self.ptcls.y(), 
                -numpy.sin(theta) * self.ptcls.x() + numpy.cos(theta) * self.ptcls.y()])

    def updateTrapConfig(self,newTrapConfig):
        self.trapConfiguration = newTrapConfig
        self.updater.trapConfiguration = newTrapConfig
        self.trapAcc.trapConfiguration = newTrapConfig

    def MakeOptimalCrystal(self, N,):
        

# vi: ts=4 sw=4

