from artiq import *
import numpy as np
import time

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
        self.attr_device("pmt_side")

        self.attr_argument("simulation_mode", BooleanValue(True))
        self.attr_argument("hist_nbins", FreeValue(60))
        self.attr_argument("hist_nmax", FreeValue(200))
        self.attr_argument("nreps", FreeValue(100))
        self.attr_argument("cool_t", FreeValue(3*ms))
        self.attr_argument("no_cool_t", FreeValue(30*ms))
        self.attr_argument("detect_t", FreeValue(1*ms))


    @kernel
    def cool_detect(self):
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
        with parallel:
            self.doppler_sw.pulse(self.detect_t)
            self.pmt_side.gate_rising(self.detect_t)
        return self.pmt_side.count()


    @kernel
    def collect_counts(self):
        for i in range(self.nreps):
            self.counts[i] = self.cool_detect()


    def run(self):
        self.counts = [0 for _ in range(self.nreps)]
        np.random.seed(int(time.time()*100)%4294967295)
        while True:
            if self.simulation_mode:
                self.counts = list(np.random.poisson(100, self.nreps))
                sleep_t = self.nreps*(2*self.cool_t + self.detect_t)
                time.sleep(sleep_t)
            else:
                self.collect_counts()
            hist, bins = np.histogram(self.counts, self.hist_nbins, range=(0, self.hist_nmax))
            hist = list(hist)
            hist_bins = list(bins)
            self.set_result("hk_default_hist_bins", hist_bins, store=False, realtime=True)
            self.set_result("hk_default_hist", hist, store=False, realtime=True)
            self.scheduler.pause()

if __name__ == "__main__":
    from artiq.frontend.artiq_run import run
    run()