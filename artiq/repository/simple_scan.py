from artiq import *
from artiq.coredevice.runtime_exceptions import RTIOUnderflow
import numpy as np
import time
import std_include

def print_underflow():
    print("Joe's RTIO underflow occured")

class SimpleScan(EnvExperiment):
    """Simple Scan

    """

    def build(self):
        self.attr_device("core")
        self.attr_device("scheduler")
        self.cool = std_include.StdCool(parent=self)
        gui_label = "Main Experiment"
        self.attr_argument("scan_time",
                           Scannable(default=LinearScan(0, 30e-3, 50)))
        self.detect = std_include.StdDetect(parent=self)

    @kernel
    def run_realtime(self):
        # self.detect.rep_start()
        # for i in range(int(self.detect.nreps)):
        #     self.core.break_realtime()
        #     self.cool.go()
        #     # experiment meat would go here if this were a fancy experiment
        #     self.detect.rep()
        # counts = self.detect.rep_done
        counts = [0]
        return counts

    def run(self):
        common_plot_x = self.set_result("common_plot_x", [], realtime=True)
        common_plot_y = self.set_result("common_plot_y", [], realtime=True)
        for t in self.scan_time:
            common_plot_x.append(t)
            counts = self.run_realtime()
            mean = sum(counts)/len(counts)
            common_plot_y.append(mean)
        self.hist_rpc(counts)

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