from artiq import *


class LEDTest(EnvExperiment):
    """Flash LED"""

    def build(self):
        self.attr_device("core")
        self.attr_device("led")

    @kernel
    def run(self):
        while True:
            delay(500*ms)
            self.led.on()
            delay(500*ms)
            self.led.off()
