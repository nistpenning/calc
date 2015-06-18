import numpy
import Ptcls
import pyopencl as cl
import pyopencl.array as cl_array
import sys
import os
import math

class ComputePotentialEnergy():

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
        src = open(absolutePathToKernels + '/calc_energy_gpu.cl', 
                'r').read()

        self.compEnergyF = cl.Program(self.ctx, src)
        try:
            self.compEnergyF.build(options = [])
        except:
              print("Error:")
              print(self.compEnergyF.get_build_info(self.ctx.devices[0],
                          cl.program_build_info.LOG))
              raise
        self.compEnergyF.compute_coulomb_acceleration.set_scalar_arg_dtypes(
                [None, None, None, None, None, 
                numpy.float32, numpy.float32, numpy.int32])

        self.compEnergyD = cl.Program(self.ctx, src)
        try:
            self.compEnergyD.build(options = [
])
        except:
              print("Error:")
              print(self.compEnergyD.get_build_info(self.ctx.devices[0],
                          cl.program_build_info.LOG))
              raise
        self.compEnergyD.compute_coulomb_acceleration.set_scalar_arg_dtypes(
                [None, None, None, None, None, 
                numpy.float64, numpy.float64, numpy.int32])

    def computeEnergy(self, x, y, z, q):

        coulombEnergy = cl_array.zero_like(q)
        xd = cl_array.to_device(self.queue, x)
        yd = cl_array.to_device(self.queue, y)
        zd = cl_array.to_device(self.queue, z)
        qd = cl_array.to_device(self.queue, q)
        prec = x.dtype
        if prec == numpy.float32:
            self.compEnergyF.calc_potential_energy(self.queue,
                    (x.size, ), None,
                    xd.data, yd.data, zd.data,
                    qd.data, coulombEnergy.data,
                    g_time_l = False)
        elif prec == numpy.float64:
            self.compEnergyD.calc_potential_energy(self.queue,
                    (x.size, ), None,
                    xd.data, yd.data, zd.data,
                    qd.data, coulombEnergy.data,
                    g_time_l = False)
        else:
            print "Unknown float type."

        return np.sum(coulombEnergy.get(self.queue))
        

        


# vi: ts=4 sw=4

