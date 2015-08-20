from artiq import *


class DDSTest(EnvExperiment):
    """idle kernel"""

    def build(self):
        self.attr_device("core")
        self.attr_device("doppler_sw")
        self.attr_device("repump_sw")
        self.attr_device("led")

    @kernel
    def run(self):
        no_cool_t = 30*ms
        cool_t = 3*ms
        while True:
            with parallel:
                self.led.on()
                self.doppler_sw.on()
                self.repump_sw.on()
            delay(cool_t)
            with parallel:
                self.led.off()
                self.doppler_sw.off()
                self.repump_sw.off()
            delay(no_cool_t)