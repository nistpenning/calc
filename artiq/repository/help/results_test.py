__author__ = 'rabi'

from artiq import *


class Auxiliary(HasEnvironment):
    def build(self):
        self.attr_device("core")

    @kernel
    def go(self):
        delay(10*ms)
        r = [1, 2, 3]
        #self.handle_results(r)
        return r

    def handle_results(self, r):
        #not possible to set results here
        # error is ValueError: Result database not present
        self.set_result("help_results_test",
            r, store=False, realtime=False)
        pass

class ResultsTest(EnvExperiment):
    """Results Test Program"""

    def build(self):
        self.attr_device("core")
        self.aux = Auxiliary(parent=self)

    @kernel
    def run(self):
        delay(1*ms)
        r = self.aux.go()
        self.handle_results(r)

    def handle_results(self, r):
        self.set_result("help_results_test",
            r, store=True, realtime=True)

    def analyze(self):
        pass
