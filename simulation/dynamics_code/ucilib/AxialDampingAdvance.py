# vi: ts=4 sw=4

import math
import numpy
import pyopencl.array as cl_array
import pyopencl as cl
import sys
import os

class AxialDampingAdvance():

    def __init__(self, ctx = None, queue = None):
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties = cl.command_queue_properties.PROFILING_ENABLE)

        absolutePathToKernels = os.path.dirname(
                os.path.realpath(__file__))
        src = open(absolutePathToKernels + '/axial_damping_advance.cl', 
                'r').read()

        self.axialDampingAdvF = cl.Program(self.ctx, src)
        try:
            self.axialDampingAdvF.build()
        except:
              print("Error:")
              print(self.axialDampingAdvF.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.axialDampingAdvF.advance_ptcls_axial_damping.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float32, numpy.int32])

        self.axialDampingAdvD = cl.Program(self.ctx, src)
        try:
            self.axialDampingAdvD.build()
        except:
              print("Error:")
              print(self.axialDampingAdvD.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.axialDampingAdvD.advance_ptcls_axial_damping.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float64, numpy.int32])

    def advancePtcls(self, xd, yd, zd, vxd, vyd, vzd, qd, md,
            dampingCoefficient, dt):
        """ 
            Dampen velocities along z to zero with rate
            dampingCoefficient. 
        """

        prec = xd.dtype
        if prec == numpy.float32:
            self.axialDampingAdvF.advance_ptcls_axial_damping(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float32(math.exp(-dampingCoefficient * dt)), 
               numpy.int32(xd.size),
               g_times_l = False)
        elif prec == numpy.float64:
            self.axialDampingAdvD.advance_ptcls_axial_damping(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float64(math.exp(-dampingCoefficient * dt)), 
               numpy.int32(xd.size),
               g_times_l = False)
        else:
            print "Unknown float type."

