import numpy
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array

class ConstantEAcc():


    def __init__(self, ctx = None, queue = None):
        self.eField = numpy.array([0, 0, 0], dtype = numpy.float32) 
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties=cl.command_queue_properties.PROFILING_ENABLE)

    def computeAcc(self, xd, yd, zd, vxd, vyd, vzd, qd, md, axd, ayd,
            azd, t, dt = None):
        """ 
           Compute acceleration due to a constant electric field.
        """
        axd += qd * self.eField[0] / md
        ayd += qd * self.eField[1] / md
        azd += qd * self.eField[2] / md

# vi: ts=4 sw=4

