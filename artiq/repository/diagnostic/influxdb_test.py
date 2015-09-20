__author__ = 'jwbritto'

import time, os

from artiq import *

class TestInfluxDB(EnvExperiment):
    """Text InfluxDB
    Measure system load and write to Results database. This shoul
    automatically be logged to InfluxDB."""

    def build(self):
        self.attr_device("scheduler")

    def run(self):
        load, _, _ = os.getloadavg()  # load avg over past 1 minute
        print(load)
        self.set_parameter("test_cpu_load", load)
        self.set_parameter("test_constant", 0.5)
        self.scheduler.submit(self.scheduler.pipeline_name,
            self.scheduler.expid,
            self.scheduler.priority, time.time() + 20, False)
