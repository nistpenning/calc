# vi: ts=4 sw=4

import math
import numpy
import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.clrandom as cl_random
import sys
import os

class CoolingAlongXAdvance():

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
        src = open(absolutePathToKernels + '/cooling_along_x_advance.cl', 
                'r').read()

        self.CoolingAlongXAdvF = cl.Program(self.ctx, src)
        try:
            self.CoolingAlongXAdvF.build()
        except:
              print("Error:")
              print(self.CoolingAlongXAdvF.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.CoolingAlongXAdvF.advance_ptcls_cooling_along_x.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None, None,
                numpy.float32, numpy.float32, numpy.float32,
                numpy.float32, numpy.float32,
                numpy.int32])

        self.CoolingAlongXAdvD = cl.Program(self.ctx, src)
        try:
            self.CoolingAlongXAdvD.build()
        except:
              print("Error:")
              print(self.CoolingAlongXAdvD.get_build_info(
                    self.ctx.devices[0],
                    cl.program_build_info.LOG))
              raise
        self.CoolingAlongXAdvD.advance_ptcls_cooling_along_x.set_scalar_arg_dtypes(
                [None, None, None, None, None, None, None, None, None,
                numpy.float64, numpy.float64, numpy.float64,
                numpy.float64, numpy.float64, 
                numpy.int32])

        self.generator = cl_random.RanluxGenerator(self.queue,
                num_work_items = 128, luxury = 1, seed = None,
                no_warmup = False, use_legacy_init = False,
                max_work_items = None)


    def advancePtcls(self, xd, yd, zd, vxd, vyd, vzd, qd, md,
            peakCoolingRate, peakDiffusionConstant, offset, width,
            dt):
        """ 
            Dampen velocities in the x direction and apply recoil kicks.
        """
        randNums = self.generator.normal(self.queue, (3 * xd.shape[0],),
                xd.dtype, sigma = 1.0)
        prec = xd.dtype
        if prec == numpy.float32:
            self.CoolingAlongXAdvD.advance_ptcls_cooling_along_x(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               randNums.data,
               qd.data, md.data,
               numpy.float32(peakCoolingRate),
               numpy.float32(peakDiffusionConstant), 
               numpy.float32(offset), 
               numpy.float32(width), 
               numpy.float32(dt), 
               numpy.int32(xd.size),
               g_times_l = False)
        elif prec == numpy.float64:
            self.CoolingAlongXAdvD.advance_ptcls_cooling_along_x(self.queue,
               (xd.size, ), None,
               xd.data, yd.data, zd.data,
               vxd.data, vyd.data, vzd.data,
               randNums.data,
               qd.data, md.data,
               numpy.float64(peakCoolingRate),
               numpy.float64(peakDiffusionConstant), 
               numpy.float64(offset), 
               numpy.float64(width), 
               numpy.float64(dt), 
               numpy.int32(xd.size),
               g_times_l = False)
        else:
            print "Unknown float type."
