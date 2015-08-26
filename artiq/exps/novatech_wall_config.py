__author__ = 'JWB'

import numpy as np

from artiq import *

global_default = {"m": 2, "gain": 0}

class NovatechWall(HasEnvironment):

    def build(self):
        self.attr_device("novatech_wall")

        cat1 = "Operate Wall"
        self.attr_argument("on_off", BooleanValue(False), cat1)
        self.attr_argument("do_scan", BooleanValue(False), cat1)
        self.attr_argument("scan", 
            Scannable(default=NoScan(325), global_min=10, global_max=1000,
                global_step=0.1, unit="kHz", ndecimals=1),
            cat1)
        self.attr_argument("scan_t", 
            NumberValue(10, unit="s", step=10, ndecimals=0, min=10, max=60), 
            cat1)

        cat2 = "Configure Wall"
        self.attr_argument("wall_order_m", 
            EnumerationValue(["1", "2", "3"]),
            cat2)
        self.attr_argument("scan_setup_t", 
            NumberValue(10, unit="s", step=10, ndecimals=0, min=0, max=60), 
            cat2)
        self.attr_argument("gain", 
            NumberValue(0.1, unit="V", step=0.1, ndecimals=2, min=0, max=0.5), 
            cat2)

    def do(self):
        pass

class NovatechAConfig(EnvExperiment):
    """Novatech Wall"""

    def build(self):
        self.wall = NovatechWall(parent=self)

    def run(self):
        pass