from artiq import *
from artiq.coredevice.runtime_exceptions import RTIOUnderflow
import numpy as np
import time
import std_include

def print_underflow():
    print("Joe's RTIO underflow occured")

class DefaultStd(EnvExperiment):
    """default experiment using std_include

    Default experiment that runs at priority 0. Cooling, histogram.
    """

    def build(self):
        self.attr_device("core")
        self.attr_device("scheduler")
        self.cool = std_include.StdCool(parent=self)
        self.attr_argument("no_cool_t",
                           NumberValue(50e-3, unit="s", ndecimals=2, step=50e-3))
        self.detect = std_include.StdDetect(parent=self)

    @kernel
    def run(self):
        while True:
            self.detect.rep_start()
            for i in range(int(self.detect.nreps)):
                self.core.break_realtime()
                self.cool.go()
                # experiment meat would go here if this were a fancy experiment
                self.detect.rep()
                delay(self.no_cool_t)
            counts = self.detect.rep_done()
            self.hist_rpc(counts)
            self.scheduler.pause()

    def hist_rpc(self, counts):
        """update histogram in GUI
        # I want this to b part of std_include.detect
        """
        hist, hist_bins = np.histogram(counts, self.detect.hist_nbins,
            range=(0, self.detect.hist_nmax))
        hist = hist.tolist()
        hist_bins = hist_bins.tolist()
        self.set_result("std_detect_hist_bins",
            hist_bins, store=False, realtime=True)
        self.set_result("std_detect_hist",
            hist, store=False, realtime=True)