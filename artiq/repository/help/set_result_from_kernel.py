from artiq import *

def print_underflow():
    print("RTIO underflow occured")

class WhatsWrongHere(EnvExperiment):
    """Should be able to set_result from @kernel decorated method. But this
    is broken with the current compiler. Here's the work around 
    until the new compiler is ready. 
    """

    def set_testval(self, val):
        self.set_result("testval", val, store=False, realtime=True)

    def build(self):
        self.attr_device("core")

    @kernel
    def run(self):
        #self.set_result("test", 0, store=False, realtime=True)
        self.set_testval(0)
