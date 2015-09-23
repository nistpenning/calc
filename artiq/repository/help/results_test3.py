from artiq import *

class ResultsTest4(EnvExperiment):
    """Results Test Program"""

    def build(self):
        self.attr_device("core")
        self.attr_argument("xs",
                   Scannable(default=LinearScan(0, 30e-3, 50)))

    @kernel
    def run_kernel(self):
        delay(1*ns)

    def run(self):
        for x in self.xs:
            self.run_kernel()
