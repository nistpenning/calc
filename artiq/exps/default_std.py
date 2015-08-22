from artiq import *
from artiq.coredevice.runtime_exceptions import RTIOUnderflow
import numpy as np
import time
from std_include import *

def print_underflow():
    print("Joe's RTIO underflow occured")

class DefaultStd(EnvExperiment, StdPrepare, StdDetect):
    """default experiment using std_include

    Default experiment that runs at priority 0. Cooling, histogram.
    """

    def build(self):
        self.attr_device("core")  
        self.attr_argument("no_cool_t", FreeValue(30*ms))
        self.attr_device("scheduler")
        pass

    @kernel
    def run(self):
        while True:
            for i in range(self.std_det_nreps):
                self.std_prepare()
                self.std_det_rep()
                delay(self.no_cool_t)
            self.std_det_rep_done()
            self.scheduler.pause()
        