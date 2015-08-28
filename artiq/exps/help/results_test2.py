__author__ = 'rabi'

from artiq import *
import numpy as np

"""Can't set numpy array as realtime result. """

class ResultsTest(EnvExperiment):
    """Results Test Program"""

    def build(self):
        self.attr_device("core")

    @kernel
    def run(self):
        delay(1*ms)

    def analyze(self):
        a = np.array([1, 1, 1])
        self.set_result("help_results_test2",
            a, store=True, realtime=True)
