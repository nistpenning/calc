__author__ = 'JWB'

import numpy as np

from artiq import *

global_default = {
"freq": 1e6,    "gain": 0, "phase": 0, "comment":""}

local_defaults = {
0: {"freq": 62.5e6, "gain": 0.5, "phase": 0, "comment": "Nallatech FGPA clock"},
3: {"freq": 3e6,    "gain": 0.5, "phase": 0, "comment": "ODF PDH"}}

class NovatechAConfig(EnvExperiment):
    """Novatech A"""

    def build(self):
        self.attr_device("novatech_a")

        for ch in range(4):
            default = local_defaults.get(ch, global_default)
            gui_group_name = "Channel {}".format(ch)

            self.comment_str = "comment_ch{}".format(ch)
            self.attr_argument(self.comment_str, 
                StringValue(default.get("comment")), 
                gui_group_name)

            self.freq_str = "freq_ch{}".format(ch)
            self.attr_argument(self.freq_str, 
                NumberValue(default.get("freq")/1e6, unit="MHz", step=1e-6, ndecimals=6, min=1e-7, max=171), 
                gui_group_name)

            self.phase_str = "phase_ch{}".format(ch)
            self.attr_argument(self.phase_str,
                NumberValue(default.get("phase"), unit="turns", step=0.25, ndecimals=2, min=0, max=1),
                gui_group_name)

            self.gain_str = "gain_ch{}".format(ch)
            self.attr_argument(self.gain_str,
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
                freq = eval("self.freq_ch{}".format(ch))
                gain = eval("self.gain_ch{}".format(ch))
                phase = eval("self.phase_ch{}".format(ch))
                self.novatech_a.set_freq(ch, freq*1e6)
                self.novatech_a.set_gain(ch, gain)
                self.novatech_a.set_phase(ch, phase)
        if self.save_state_to_eeprom:
            self.novatech_a.save_state_to_eeprom()

