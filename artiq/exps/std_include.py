from artiq import *
import numpy as np

class StdPrepare():
    """Standard qubit state perparation
    """
    def __init__(self):
        super(StdPrepare, self).__init__()
    def build(self):
        super(StdPrepare, self).build()  
        self.attr_device("doppler_sw")
        self.attr_device("repump_sw")
        self.attr_argument("std_prepare_cool_t", FreeValue(3*ms))
        self.attr_argument("std_prepare_repump_t", FreeValue(1*ms))
          
    @kernel
    def std_prepare(self):
        """prepare qubit by cooling and repumping
        """
        with parallel:
            self.led.on()
            self.doppler_sw.on()
            self.repump_sw.on() 
            delay(self.std_prepare_cool_t)
            self.doppler_sw.off()
            delay(self.std_prepare_repump_t) 
            self.repump_sw.off()
            self.led.off()      

class StdDetect():
    """Standard qubit state detection
    """
    def __init__(self):
        super(StdDetect, self).__init__()
    def build(self):
        super(StdDetect, self).build()
        self.attr_device("doppler_sw")
        self.attr_device("repump_sw")
        self.attr_device("led")
        self.attr_device("sideview_pmt")

        self.attr_argument("std_det_hist_nbins", FreeValue(60))
        self.attr_argument("std_det_hist_nmax", FreeValue(200))
        self.attr_argument("std_det_nreps", FreeValue(100))
        self.attr_argument("std_detect_t", FreeValue(1*ms))

        self.std_detect_counts = [0 for _ in range(self.std_det_nreps)]
        self.std_detect_rep_i = 0

    @kernel
    def std_detect_rep(self):
        """projective measurement; of single experiment repetition
        """
        self.core.break_realtime()
        with parallel:
            self.led.pulse(self.std_detect_t)
            self.doppler_sw.pulse(self.std_detect_t)
            self.sideview_pmt.gate_rising(self.std_detect_t)
        cts = self.sideview_pmt.count()
        self.std_detect_counts[self.std_detect_rep_i] = cts
        self.std_detect_rep_i += 1

    @kernel
    def std_detect_rep_done(self):
        """call at conclusion of a set of identical experiment reps
        """
        self.std_detect_counts = [0 for _ in range(self.std_det_nreps)]
        self.std_detect_rep_i = 0
        mean, std = self.std_detect_rep_hist_rpc()
        return mean, std

    def std_detect_rep_hist_rpc(self):
        """update histogram in GUI
        """
        counts = self.std_detect_counts
        hist, bins = np.histogram(counts, self.std_det_hist_nbins, 
            range=(0, self.std_det_hist_nmax))
        hist = list(hist)
        hist_bins = list(bins)
        self.set_result("std_det_hist_bins", 
            hist_bins, store=False, realtime=True)
        self.set_result("std_det_hist", 
            hist, store=False, realtime=True)
        mean = np.mean(counts)
        std = np.std(counts)
        return mean, std