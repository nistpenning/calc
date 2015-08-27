__author__ = 'JWB'

import numpy as np

from artiq import *

global_default = {
"freq": 1e6,    "gain": 0, "phase": 0, "comment":""}

local_defaults = {
0: {"freq": 62.5e6, "gain": 0.5, "phase": 0, "comment": "Nallatech FGPA clock"},
1: {"freq": 125e6, "gain": 0.5, "phase": 0, "comment": "KC705 clock"},
3: {"freq": 10e6,    "gain": 0.5, "phase": 0, "comment": "ODF PDH"}}

class NovatechAConfig(EnvExperiment):
    """Novatech A"""

    def build(self):
        self.attr_device("novatech_a")

        for ch in range(4):
            default = local_defaults.get(ch, global_default)
            gui_group_name = "Channel {}".format(ch)

            self.get_argument("comment_ch{}".format(ch), 
                StringValue(default.get("comment")), 
                gui_group_name)

            self.get_argument("freq_ch{}".format(ch), 
                NumberValue(default.get("freq")/1e6, unit="MHz", step=1e-6, ndecimals=6, min=1e-7, max=171), 
                gui_group_name)

            self.get_argument("phase_ch{}".format(ch),
                NumberValue(default.get("phase"), unit="turns", step=0.25, ndecimals=2, min=0, max=1),
                gui_group_name)

            self.get_argument("gain_ch{}".format(ch),
                NumberValue(default.get("gain"), unit="Volts", step=0.1, ndecimals=1, min=0, max=0.50),
                gui_group_name)  

        gn2 = "Configure Action(s) Upon Run"
        self.attr_argument("reset", BooleanValue(False), gn2)
        self.attr_argument("set_freq_phase_gain", BooleanValue(False), gn2)
        self.attr_argument("save_state_to_eeprom", BooleanValue(False), gn2)

    def run(self):
        if self.reset:
            self.novatech_a.reset()
        if self.set_freq_phase_gain:
            for ch in range(4):
                freq = getattr(self, "freq_ch{}".format(ch))
                gain = getattr(self,  "gain_ch{}".format(ch))
                phase = getattr(self, "phase_ch{}".format(ch))
                self.novatech_a.set_freq(ch, freq*1e6)
                self.novatech_a.set_gain(ch, gain)
                self.novatech_a.set_phase(ch, phase)
        if self.save_state_to_eeprom:
            self.novatech_a.save_state_to_eeprom()

