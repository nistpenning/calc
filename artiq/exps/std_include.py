from artiq import *
import numpy as np

class StdCool(HasEnvironment):
    """Standard qubit cooling and state perparation
    """
    def build(self):
        self.attr_device("core")
        self.attr_device("doppler_sw")
        self.attr_device("repump_sw")
        self.attr_device("led")

        gui_group = "Standard Cool"
        self.attr_argument("cool_t",
                           NumberValue(3e-3, unit="s", step=1e-3, ndecimals=3, min=0, max=1),
                           gui_group)
        self.attr_argument("repump_t",
                           NumberValue(1e-3, unit="s", step=1e-3, ndecimals=3, min=0, max=1),
                           gui_group)
          
    @kernel
    def go(self):
        """prepare qubit by cooling and repumping
        """
        with parallel:
            self.led.on()
            self.doppler_sw.on()
            self.repump_sw.on() 
            delay(self.cool_t)
            self.doppler_sw.off()
            delay(self.repump_t)
            self.repump_sw.off()
            self.led.off()

class StdDetect(HasEnvironment):
    """Standard qubit state detection
    """
    def build(self):
        self.attr_device("core")
        self.attr_device("doppler_sw")
        self.attr_device("repump_sw")
        self.attr_device("led")
        self.attr_device("sideview_pmt")

        gui_group = "Standard Detect"
        self.attr_argument("hist_nbins",
                           NumberValue(60, step=10, ndecimals=0, min=0),
                           gui_group)
        self.attr_argument("hist_nmax",
                           NumberValue(200, step=10, ndecimals=0, min=10),
                           gui_group)
        self.attr_argument("nreps",
                           NumberValue(100, step=50, ndecimals=0, min=1, max=1000),
                           gui_group)
        self.attr_argument("detect_t",
                           NumberValue(1e-3, unit="s", step=1e-3, ndecimals=4, min=1e-6, max=1),
                           gui_group)

    @kernel
    def rep(self):
        """projective measurement; of single experiment repetition
        """
        with parallel:
            self.led.pulse(self.detect_t)
            self.doppler_sw.pulse(self.detect_t)
            self.sideview_pmt.gate_rising(self.detect_t)
        cts = self.sideview_pmt.count()
        self.counts[self.rep_i] = cts
        self.rep_i += 1

    @kernel
    def rep_done(self):
        """call at conclusion of a set of identical experiment reps
        """
        return self.counts

    @kernel
    def rep_start(self):
        """initialize count buffer
        """
        self.rep_i = 0
        self.counts = [0 for _ in range(int(self.nreps))]

