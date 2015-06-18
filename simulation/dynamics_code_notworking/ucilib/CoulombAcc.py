import numpy
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array
import sys
import os
import math


class CoulombAcc():

    BLOCK_SIZE = 256
    PTCL_UNROLL_FACTOR = 1

    maxNumThreadsX = 2**12

    k = 1./(4.*numpy.pi*8.854187817620e-12)
    impactFact = 1.0e-9**2 

    def __init__(self, ctx = None, queue = None):
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.mf = cl.mem_flags

        absolutePathToKernels = os.path.dirname(
                os.path.realpath(__file__))
        src = open(absolutePathToKernels + '/compute_coulomb_acceleration.cl', 
                'r').read()

        self.compAccF = cl.Program(self.ctx, src)
        try:
            self.compAccF.build(options = [
              ' -DBLOCK_SIZE=' + str(self.BLOCK_SIZE) +
              ' -DPTCL_UNROLL_FACTOR=' + str(self.PTCL_UNROLL_FACTOR)])
        except:
              print("Error:")
              print(self.compAccF.get_build_info(self.ctx.devices[0],
                          cl.program_build_info.LOG))
              raise
        self.compAccF.compute_coulomb_acceleration.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float32, numpy.float32, numpy.int32, None, None,
                None])

        self.compAccD = cl.Program(self.ctx, src)
        try:
            self.compAccD.build(options = [
              ' -DUSE_DOUBLE=TRUE -DBLOCK_SIZE=' + str(self.BLOCK_SIZE) +
              ' -DPTCL_UNROLL_FACTOR=' + str(self.PTCL_UNROLL_FACTOR)])
        except:
              print("Error:")
              print(self.compAccD.get_build_info(self.ctx.devices[0],
                          cl.program_build_info.LOG))
              raise
        self.compAccD.compute_coulomb_acceleration.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float64, numpy.float64, numpy.int32, None, None,
                None])


    def computeAcc(self, xd, yd, zd, vxd, vyd, vzd, qd, md, axd, ayd,
            azd, t, dt = None):
        """ 
           Compute acceleration due to coulomb forces between particles.
           All variables are device arrays.  They are assumed to be the
           same length and precision.  The computed accelerations are
           accumulated into the a?d arrays (i.e. the a?d arrays are not
           set to zero).
        """

        self.computeLaunchConfig(xd.size)

        prec = xd.dtype
        if prec == numpy.float32:
            self.compAccF.compute_coulomb_acceleration(self.queue,
               (self.numThreadsX, self.numThreadsY),
               (self.BLOCK_SIZE, 1),
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float32(self.k), numpy.float32(self.impactFact), 
               numpy.int32(xd.size),
               axd.data, ayd.data, azd.data,
               g_times_l = False)
        elif prec == numpy.float64:
            self.compAccD.compute_coulomb_acceleration(self.queue,
               (self.numThreadsX, self.numThreadsY),
               (self.BLOCK_SIZE, 1),
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float64(self.k), numpy.float64(self.impactFact), 
               numpy.int32(xd.size),
               axd.data, ayd.data, azd.data,
               g_times_l = False)
        else:
            print "Unknown float type."


    def computeLaunchConfig(self, numPtcls):
        self.numThreads = int(math.ceil(float(numPtcls) / self.PTCL_UNROLL_FACTOR))
        self.numThreadsX = self.BLOCK_SIZE * int(math.ceil(float(self.numThreads) / 
                        float(self.BLOCK_SIZE)))
        if self.numThreadsX > self.maxNumThreadsX:
            self.numThreadsX = self.maxNumThreadsX
        self.numThreadsY = int(math.ceil(float(self.numThreads) / 
                    float(self.numThreadsX)))


# vi: ts=4 sw=4

