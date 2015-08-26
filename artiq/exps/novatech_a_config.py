__author__ = 'JWB'

import numpy as np

from artiq import *

global_default = {
"freq": 1e6,    "gain": 0, "phase": 0, "comment":""}

local_defaults = {
0: {"freq": 62.5e6, "gain": 0.5, "phase": 0, "comment": "Nallatech FGPA clock"},
3: {"freq": 3e6,    "gain": 0.5, "phase": 0, "comment": "ODF PDH"}}

class NovatechChannel(HasEnvironment):
    def __init__(self, parent, ch):
        self.parent = parent
        self.ch = ch
        super().__init__(parent=parent)

    def build(self):
        self.attr_device("novatech_a")

        default = local_defaults.get(self.ch, global_default)
        gui_group_name = "Channel {}".format(self.ch)

        self.comment_str = "comment_ch{}".format(self.ch)
        self.attr_argument(self.comment_str, 
            StringValue(default.get("comment")), 
            gui_group_name)

        self.freq_str = "freq_ch{}".format(self.ch)
        self.attr_argument(self.freq_str, 
            NumberValue(default.get("freq")/1e6, unit="MHz", step=1e-6, ndecimals=6, min=1e-7, max=171), 
            gui_group_name)

        self.gain_str = "gain_ch{}".format(self.ch)
        self.attr_argument(self.gain_str,
            NumberValue(default.get("gain"), unit="Volts", step=0.1, ndecimals=1, min=0, max=0.50),
            gui_group_name)  

        self.phase_str = "phase_ch{}".format(self.ch)
        self.attr_argument(self.phase_str,
            NumberValue(default.get("phase"), unit="turns", step=0.25, ndecimals=2, min=0, max=1),
            gui_group_name)

    def do(self):
        freq = eval("self.{}".format(self.freq_str))
        gain = eval("self.{}".format(self.gain_str))
        phase = eval("self.{}".format(self.phase_str))
        self.novatech_a.set_freq(self.ch, freq*1e6)
        self.novatech_a.set_gain(self.ch, gain)
        self.novatech_a.set_phase(self.ch, phase)

class NovatechAConfig(EnvExperiment):
    """Novatech A"""

    def build(self):
        self.ch0 = NovatechChannel(parent=self, ch=0)
        self.ch1 = NovatechChannel(parent=self, ch=1)
        self.ch2 = NovatechChannel(parent=self, ch=2)
        self.ch3 = NovatechChannel(parent=self, ch=3)

    def run(self):
        self.ch0.do()
        self.ch1.do()
        self.ch2.do()
        self.ch3.do()