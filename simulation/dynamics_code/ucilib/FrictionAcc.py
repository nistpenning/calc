# vi: ts=4 sw=4

import numpy
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array
import sys
import os
import math


class FrictionAcc():

# Friction coefficient in N/(m/s)
    def __init__(self, ctx = None, queue = None):
        self.frictionCoefficient = 0

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
        src = open(absolutePathToKernels + '/compute_fric_acc.cl', 
                'r').read()

        self.compAccF = cl.Program(self.ctx, src)
        try:
            self.compAccF.build(options = [])
        except:
              print("Error:")
              print(self.compAccF.get_build_info(self.ctx.devices[0],
                          cl.program_build_info.LOG))
              raise

        self.compAccD = cl.Program(self.ctx, src)
        try:
            self.compAccD.build(options = [' -DUSE_DOUBLE=TRUE'])
        except:
              print("Error:")
              print(self.compAccD.get_build_info(self.ctx.devices[0],
                          cl.program_build_info.LOG))
              raise


    def computeAcc(self, xd, yd, zd, vxd, vyd, vzd, qd, md, axd, ayd,
            azd, t, dt = None):
        """ 
           Compute acceleration due to friction.
        """
        
        prec = xd.dtype
        if prec == numpy.float32:
            self.compAccF.compute_acceleration(self.queue,
               (xd.size, ),
               None,
               vxd.data, vyd.data, vzd.data, md.data,
               numpy.float32(self.frictionCoefficient),
               numpy.int32(xd.size),
               axd.data, ayd.data, azd.data,
               g_times_l = False)
        elif prec == numpy.float64:
            self.compAccF.compute_acceleration(self.queue,
               (xd.size, ),
               None,
               vxd.data, vyd.data, vzd.data,
               md.data,
               numpy.float64(self.frictionCoefficient),
               numpy.int32(xd.size),
               axd.data, ayd.data, azd.data,
               g_times_l = False)
        else :
            print "Unknown float type."


