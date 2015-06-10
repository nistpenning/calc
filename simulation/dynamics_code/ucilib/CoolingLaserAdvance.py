
# vi: ts=4 sw=4

import math
import numpy
import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.clrandom as cl_random
import sys
import os

class CoolingLaserAdvance():


    def __init__(self, ctx = None, queue = None):#, displacement):

        self._PlanckConstantReduced = 1.0545717e-34
        # wavelength of cooling laser
        lam = 313.0e-9
        # wave vector
        self.k0 = numpy.array([0, 0, 2.0 * numpy.pi / lam],
                dtype = numpy.float32)
        self.x0 = numpy.array([0, 0, 0], dtype = numpy.float32)
        # 1/e radius of cooling laser
        self.sigma = 1.0e-3

        # line width (unsaturated)
        self.gamma = 2.0 * numpy.pi * 19.0e6
        # Detuning at zero velocity
        self.delta0 = -0.5 * self.gamma
        # Saturation parameter
        self.S = 0.1
        # Offset from origin
        #self.d = numpy.array([0, displacement, 0], dtype = numpy.float32)



        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties = cl.command_queue_properties.PROFILING_ENABLE)
        absolutePathToKernels = os.path.dirname(
                os.path.realpath(__file__))
        src = open(absolutePathToKernels +
                '/cooling_laser_advance.cl', 'r').read()
        self.program = cl.Program(self.ctx, src)
        try:
            self.program.build()
        except:
              print("Error:")
              print(self.program.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.program.compute_mean_scattered_photons_homogeneous_beam.set_scalar_arg_dtypes(
                [None, None, None, None, None, None,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.float32,
                numpy.int32,
                None])
        self.program.compute_mean_scattered_photons_gaussian_beam.set_scalar_arg_dtypes(
                [None, None, None, None, None, None,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.float32, numpy.float32, numpy.float32, numpy.float32,
                numpy.float32,
                numpy.int32,
                None])
        self.program.countEmissions.set_scalar_arg_dtypes(
                [None, None, numpy.int32, None, numpy.int32]
                )
        self.program.computeKicks.set_scalar_arg_dtypes(
                [None, None, numpy.int32, None,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.float32, 
                numpy.float32,
                None, None, None,
                numpy.int32])

        self.generator = cl_random.RanluxGenerator(self.queue,
                num_work_items = 128, luxury = 1, seed = None,
                no_warmup = False, use_legacy_init = False,
                max_work_items = None)


    def computeAcc(self, xd, yd, zd, vxd, vyd, vzd, qd, md, 
            axd, ayd, azd, t, dt):

        # Compute average numbers of scattered photons
        nbars = cl_array.zeros_like(xd)
        if self.sigma == None:
            self.program.compute_mean_scattered_photons_homogeneous_beam(
                    self.queue, (xd.size, ), None,
                    xd.data, yd.data, zd.data,
                    vxd.data, vyd.data, vzd.data,
                    numpy.float32(self.k0[0]),
                    numpy.float32(self.k0[1]),
                    numpy.float32(self.k0[2]),
                    numpy.float32(self.gamma),
                    numpy.float32(self.delta0),
                    numpy.float32(self.S),
                    numpy.float32(dt),
                    numpy.int32(xd.size),
                    nbars.data)
        else:
            self.program.compute_mean_scattered_photons_gaussian_beam(
                    self.queue, (xd.size, ), None,
                    xd.data, yd.data, zd.data,
                    vxd.data, vyd.data, vzd.data,
                    numpy.float32(self.k0[0]),
                    numpy.float32(self.k0[1]),
                    numpy.float32(self.k0[2]),
                    numpy.float32(self.x0[0]),
                    numpy.float32(self.x0[1]),
                    numpy.float32(self.x0[2]),
                    numpy.float32(self.sigma),
                    numpy.float32(self.gamma),
                    numpy.float32(self.delta0),
                    numpy.float32(self.S),
                    numpy.float32(dt),
                    numpy.int32(xd.size),
                    nbars.data)
        
        # Compute scattered photons and associated recoil kicks
        nMax = int(math.ceil(10.0 * self.S * 
                    (self.gamma / 2.0 / numpy.pi) * dt))
#        print nMax
#        print (10.0 * self.S * 
#                    (self.gamma / 2.0 / numpy.pi) * dt)
        actualNs = self.findSample(nbars, nMax)
        recoilDirectionsD = cl_array.Array(self.queue,
                [nbars.size, nMax, 3], dtype = numpy.float32)
        self.generator.fill_normal(recoilDirectionsD)

        # apply recoil kicks to particles
        recoilMomentum = numpy.linalg.norm(self.k0) * self._PlanckConstantReduced
        self.program.computeKicks(
                    self.queue, (xd.size, ), None,
                    md.data,
                    actualNs.data,
                    numpy.int32(nMax),
                    recoilDirectionsD.data,
                    numpy.float32(self.k0[0]),
                    numpy.float32(self.k0[1]),
                    numpy.float32(self.k0[2]),
                    numpy.float32(recoilMomentum),
                    numpy.float32(dt),
                    axd.data, ayd.data, azd.data,
                    numpy.int32(xd.shape[0]))


    def findSample(self, meanN, nMax):
        auxRandNumbers = cl_array.Array(self.queue,
                [meanN.shape[0], nMax], numpy.float32)
        self.generator.fill_uniform(auxRandNumbers)
        sample = cl_array.Array(self.queue, [meanN.shape[0],],
                dtype = numpy.int32)
        self.program.countEmissions(
                self.queue, (meanN.size, ), None,
                meanN.data, auxRandNumbers.data, numpy.int32(nMax),
                sample.data, numpy.int32(meanN.shape[0]))
        return sample


