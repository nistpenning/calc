from artiq import *

class CountsTest(EnvExperiment):
    """Counts Test"""

    def build(self):
            self.attr_device("core")
            self.attr_device("pmt0")
            self.attr_device("ttl1")  # doppler_sw
            self.attr_device("ttl2")  # repump_sw
            self.attr_device("ttl12")  # sideview PMT/Andor switch
    @kernel
    def run(self):
        while True:
            self.core.break_realtime()
            delay(1*s)
            self.ttl12.on()
            self.ttl1.on()
            self.ttl2.on()
            self.pmt0.gate_rising(1000*us)
            delay(1100*us)
            print(self.pmt0.count())