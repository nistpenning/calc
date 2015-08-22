from artiq import *
from artiq.coredevice.runtime_exceptions import RTIOUnderflow
import numpy as np
import time

def print_underflow():
    print("Joe's RTIO underflow occured")

class Default(EnvExperiment):
    """default

    Default experiment that runs at priority 0. Cooling, histogram.
    """

    def build(self):
        self.attr_device("core")
        self.attr_device("scheduler")
        self.attr_device("doppler_sw")
        self.attr_device("repump_sw")
        self.attr_device("led")
        self.attr_device("sideview_pmt")
        self.attr_argument("no_cool_t", FreeValue(30*ms))

        self.attr_argument("simulation_mode", BooleanValue(False))
        self.attr_argument("hist_nbins", FreeValue(60))
        self.attr_argument("hist_nmax", FreeValue(200))
        self.attr_argument("nreps", FreeValue(100))
        self.attr_argument("cool_t", FreeValue(3*ms))
        self.attr_argument("detect_t", FreeValue(1*ms))

    @kernel
    def cool_detect(self):
        self.core.break_realtime()
        with parallel:
            self.led.on()
            self.doppler_sw.on()
            self.repump_sw.on()
        delay(self.cool_t)
        with parallel:
            self.led.off()
            self.doppler_sw.off()
            self.repump_sw.off()
        delay(self.no_cool_t)
        try: 
            with parallel:
                self.doppler_sw.pulse(self.detect_t)
                self.sideview_pmt.gate_rising(self.detect_t)
            a = self.sideview_pmt.count()
        except RTIOUnderflow:
            print_underflow()
        return a

    @kernel
    def run(self):
        while True:
            counts = [0 for _ in range(self.nreps)]
            if self.simulation_mode:
                counts = list(np.random.poisson(100, self.nreps))
                sleep_t = self.nreps*(2*self.cool_t + self.detect_t)
                time.sleep(sleep_t)
            else:
                for i in range(self.nreps):
                    counts[i] = self.cool_detect()
            self.build_hist(counts)
            delay(self.no_cool_t)
        
    def build_hist(self, counts):
        # build the histogram on the master
        hist, bins = np.histogram(counts, self.hist_nbins, 
            range=(0, self.hist_nmax))
        hist = list(hist)
        hist_bins = list(bins)
        self.set_result("hk_default_hist_bins", 
            hist_bins, store=False, realtime=True)
        self.set_result("hk_default_hist", 
            hist, store=False, realtime=True)

if __name__ == "__main__":
    from artiq.frontend.artiq_run import run
    run()