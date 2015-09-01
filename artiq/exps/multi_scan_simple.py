from artiq import *
from artiq.coredevice.runtime_exceptions import RTIOUnderflow
import numpy as np
import time

def print_underflow():
    print("Joe's RTIO underflow occured")

class MultiScan(EnvExperiment):
    """multi scan

    Default experiment that runs at priority 0. Cooling, histogram.
    """

    def build(self):
        self.attr_argument("x1", NumberValue(3))
        self.attr_argument("x2", NumberValue(3))
        self.attr_argument("xs",
                   Scannable(default=LinearScan(0, 1, 5)))
        scan_param_options = ["x1", "x2"]
        self.attr_argument("scan_param",
            EnumerationValue(scan_param_options, "x1"))


    def multif(self, x1, x2):
        print(x1,x2)

    def run(self):
        for xi in self.xs:
            setattr(self, self.scan_param, xi)
            self.multif(self.x1, self.x2)

