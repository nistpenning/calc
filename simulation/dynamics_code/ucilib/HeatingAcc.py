import numpy
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array
import pyopencl.clrandom as cl_random
import sys
import os

class HeatingAcc():

    def __init__(self, ctx = None, queue = None):
        self.diffusionConstant = 0.0
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties=cl.command_queue_properties.PROFILING_ENABLE)
        self.mf = cl.mem_flags

        self.generator = cl_random.RanluxGenerator(self.queue,
                num_work_items = 128, luxury = 1, seed = None,
                no_warmup = False, use_legacy_init = False,
                max_work_items = None)

    def computeAcc(self, xd, yd, zd, vxd, vyd, vzd, qd, md, axd, ayd,
            azd, t, dt):
#        print "sigma: ", (1.0 / 3.0) * self.diffusionConstant / numpy.sqrt(dt)
        axd += self.generator.normal(self.queue, axd.shape, axd.dtype,
                sigma = (1.0 / 3.0) * self.diffusionConstant / numpy.sqrt(dt));
        ayd += self.generator.normal(self.queue, ayd.shape, ayd.dtype,
                sigma = (1.0 / 3.0) * self.diffusionConstant / numpy.sqrt(dt));
        azd += self.generator.normal(self.queue, azd.shape, azd.dtype,
                sigma = (1.0 / 3.0) * self.diffusionConstant / numpy.sqrt(dt));

# vi: ts=4 sw=4

