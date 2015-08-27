__author__ = 'JWB'

import numpy as np

from artiq import *

global_default = {
"freq": 1e6,    "gain": 0, "phase": 0, "comment":""}

local_defaults = {
0: {"freq": 62.5e6, "gain": 0.5, "phase": 0, "comment": "comment"},
1: {"freq": 125e6, "gain": 0.5, "phase": 0, "comment": ""},
3: {"freq": 3e6,    "gain": 0.5, "phase": 0, "comment": ""}}

class CoreDDSDebug(EnvExperiment):
    """Core DDS: Set & Debug"""
    numdds = 8

    def build(self):
        self.attr_device("core")
        [self.attr_device("dds{}".format(x)) for x in range(self.numdds)]

        def make_gui(ch):
            default = local_defaults.get(ch, global_default)
            gui_group_name = "Channel {}".format(ch)

            self.attr_argument("freq_ch{}".format(ch), 
                NumberValue(default.get("freq")/1e6, unit="MHz", step=1e-6, ndecimals=6, min=1e-7, max=171), 
                gui_group_name)

            self.attr_argument("phase_ch{}".format(ch),
                NumberValue(default.get("phase"), unit="turns", step=0.25, ndecimals=2, min=0, max=1),
                gui_group_name)

        for ch in range(self.numdds):
            make_gui(ch)
        make_gui("ALL")

        gn2 = "Configure Action(s) Upon Run"
        chnos = ["{}".format(x) for x in range(self.numdds)]
        chnos.append("ALL")
        chnos.append("NONE")
        self.attr_argument("set_dds", 
            EnumerationValue(chnos, "ALL"), 
            gn2)
        self.attr_argument("reset", 
            EnumerationValue(chnos, "NONE"), 
            gn2)

    @kernel
    def run_kernel(self, ch, freq, phase):
        self.core.break_realtime()
        fn = "dds{}".format(ch)
        f = getattr(self, fn)
        f.set(freq*MHz, phase)
    
    def run(self):
        if self.reset == "NONE":
            chnos = []
        elif self.reset == "ALL":
            chnos = range(self.numdds)
        else:
            chnos = [int(self.reset)]
        for ch in chnos:
            # reset DDS (not yet implemented)
            pass

        if self.set_dds == "NONE":
            chnos = []
        elif self.set_dds == "ALL":
            chnos = range(self.numdds)
        else:
            chnos = [int(self.set_dds)]
        print("chnos = {}".format(chnos))
        for ch in chnos:
            if self.set_dds == "ALL":
                freq = getattr(self, "freq_chALL")
                phase = getattr(self, "phase_chALL")
            else:
                freq = getattr(self, "freq_ch{}".format(ch))
                phase = getattr(self, "phase_ch{}".format(ch))
            self.run_kernel(ch, freq, phase)
