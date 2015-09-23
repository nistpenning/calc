from artiq import *

class TestTTLout(EnvExperiment):
    """Test all TTL outputs"""
    def build(self):
        self.attr_device("core")
        self.attr_device("dds_bus")
        self.attr_device("ttl0")
        self.attr_device("ttl1")
        self.attr_device("ttl2")
        self.attr_device("ttl3")
        self.attr_device("ttl4")
        self.attr_device("ttl5")
        self.attr_device("ttl6")
        self.attr_device("ttl7")
        self.attr_device("ttl8")
        self.attr_device("ttl9")
        self.attr_device("ttl10")
        self.attr_device("ttl11")
        self.attr_device("ttl12")
        self.attr_device("ttl13")
        self.attr_device("ttl14")
        self.attr_device("ttl15")
        self.pt = 5*us

    @kernel
    def run(self):
        while True:
            self.ttl0.pulse(self.pt)
            self.ttl1.pulse(self.pt)
            self.ttl2.pulse(self.pt)
            self.ttl3.pulse(self.pt)
            self.ttl4.pulse(self.pt)
            self.ttl5.pulse(self.pt)
            self.ttl6.pulse(self.pt)
            self.ttl7.pulse(self.pt)
            self.ttl8.pulse(self.pt)
            self.ttl9.pulse(self.pt)
            self.ttl10.pulse(self.pt)
            self.ttl11.pulse(self.pt)
            self.ttl12.pulse(self.pt)
            self.ttl13.pulse(self.pt)
            self.ttl14.pulse(self.pt)
            self.ttl15.pulse(self.pt)
            delay(20*ms)
