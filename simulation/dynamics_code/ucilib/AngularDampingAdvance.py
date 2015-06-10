# vi: ts=4 sw=4

import math
import numpy
import pyopencl.array as cl_array
import pyopencl as cl
import sys
import os

class AngularDampingAdvance():

    def __init__(self, ctx = None, queue = None):
        self.minRadius = 1.0e-5
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties = cl.command_queue_properties.PROFILING_ENABLE)
        absolutePathToKernels = os.path.dirname(
                os.path.realpath(__file__))
        src = open(absolutePathToKernels + '/angular_damping_advance.cl', 
                'r').read()

        self.angularDampingAdvF = cl.Program(self.ctx, src)
        try:
            self.angularDampingAdvF.build()
        except:
              print("Error:")
              print(self.angularDampingAdvF.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.angularDampingAdvF.advance_ptcls_angular_damping.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.int32])

        self.angularDampingAdvD = cl.Program(self.ctx, src)
        try:
            self.angularDampingAdvD.build()
        except:
              print("Error:")
              print(self.angularDampingAdvD.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.angularDampingAdvD.advance_ptcls_angular_damping.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float64, numpy.float64, numpy.float64,
                numpy.int32])

    def advancePtcls(self, xd, yd, zd, vxd, vyd, vzd, qd, md,
            dampingCoefficient, omega, dt):
        """ 
            Dampen velocities in the x-y plane.
        """
        prec = xd.dtype
        if prec == numpy.float32:
            self.angularDampingAdvD.advance_ptcls_angular_damping(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float32(math.exp(-dampingCoefficient * dt)),
               numpy.float32(omega), 
               numpy.float32(self.minRadius), 
               numpy.int32(xd.size),
               g_times_l = False)
        elif prec == numpy.float64:
            self.angularDampingAdvD.advance_ptcls_angular_damping(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float64(math.exp(-dampingCoefficient * dt)),
               numpy.float64(omega), 
               numpy.float64(self.minRadius), 
               numpy.int32(xd.size),
               g_times_l = False)
        else:
            print "Unknown float type."
