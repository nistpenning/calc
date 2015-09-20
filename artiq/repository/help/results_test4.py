from artiq import *

class ResultsTest4(EnvExperiment):
    """Results Test Program"""

    def build(self):
        self.attr_device("core")
        self.attr_argument("xs",
                   Scannable(default=LinearScan(0, 30e-3, 100)))

    @kernel
    def run_kernel(self):
        delay(1*ns)
        dd = mu_to_seconds(now_mu())
        return t

    def run(self):
        ts = []
        for x in self.xs:
            t = self.run_kernel()
            ts.append(t)
        self.handle_result(ts)

    def handle_result(self, ts):
        self.set_result("ts", ts)
