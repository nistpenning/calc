import numpy
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array
import sys
import os
import math


class CyclAdvance():

    def __init__(self, ctx = None, queue = None):
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties=cl.command_queue_properties.PROFILING_ENABLE)

        absolutePathToKernels = os.path.dirname(
                os.path.realpath(__file__))
        src = open(absolutePathToKernels + '/cycl_advance.cl', 
                'r').read()

        self.cyclAdvF = cl.Program(self.ctx, src)
        try:
            self.cyclAdvF.build()
        except:
              print("Error:")
              print(self.cyclAdvF.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.cyclAdvF.advance_ptcls_cycl.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float32, numpy.float32, numpy.int32])

        self.cyclAdvD = cl.Program(self.ctx, src)
        try:
            self.cyclAdvD.build()
        except:
              print("Error:")
              print(self.cyclAdvD.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.cyclAdvD.advance_ptcls_cycl.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None,
                numpy.float64, numpy.float64, numpy.int32])


    def advancePtcls(self, xd, yd, zd, vxd, vyd, vzd, qd, md, Bz, dt):
        """ 
           Advance particles for time dt in a uniform magnetic field of
           magnitude Bz along the z direction.  The z position advanced
           according to vzd.
        """

        prec = xd.dtype
        if prec == numpy.float32:
            self.cyclAdvF.advance_ptcls_cycl(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float32(Bz), numpy.float32(dt), 
               numpy.int32(xd.size),
               g_times_l = False)
        elif prec == numpy.float64:
            self.cyclAdvD.advance_ptcls_cycl(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               qd.data, md.data,
               numpy.float64(Bz), numpy.float64(dt), 
               numpy.int32(xd.size),
               g_times_l = False)
        else:
            print "Unknown float type."


# vi: ts=4 sw=4

