__author__ = 'jwbritto'
from artiq import *
import time


class PreemptTest(EnvExperiment):
    """Preempt test"""
    def build(self):
        self.attr_device("scheduler")
        pass

    def run(self):
        while True:
            time.sleep(2)
            self.scheduler.pause()