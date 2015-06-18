# vi: ts=4 sw=4

import numpy
import ptcls
import pyopencl as cl
import pyopencl.array as cl_array
import sys
import os


class BorisUpdater():

    def __init__(self, ctx = None, queue = None):
        self.ctx = ctx
        self.queue = queue
        if self.ctx == None:
            self.ctx = cl.create_some_context()
        if self.queue == None:
            self.queue = cl.CommandQueue(self.ctx,
                properties=cl.command_queue_properties.PROFILING_ENABLE)

    def update(self, xd, yd, zd, vxd, vyd, vzd, qd, md, accelerations, t, dt):

# First half kick
        axd = cl_array.zeros_like(xd)
        ayd = cl_array.zeros_like(xd)
        azd = cl_array.zeros_like(xd)
        for acc in accelerations:
            acc.computeAcc(xd, yd, zd, vxd, vyd, vzd, qd, md,
                    axd, ayd, azd, t)
        vxd += (0.5 * dt) * axd
        vyd += (0.5 * dt) * ayd
        vzd += (0.5 * dt) * azd
        

# Advance positions
        xd += dt * vxd
        yd += dt * vyd
        zd += dt * vzd
        
#Second half kick
        
        axd.fill(0.0, self.queue)
        ayd.fill(0.0, self.queue)
        azd.fill(0.0, self.queue)
        for acc in accelerations:
            acc.computeAcc(xd, yd, zd, vxd, vyd, vzd, qd, md,
                    axd, ayd, azd, t)
        vxd += (0.5 * dt) * axd
        vyd += (0.5 * dt) * ayd
        vzd += (0.5 * dt) * azd

