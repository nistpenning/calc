from artiq import *
import logging

from artiq.tools import verbosity_args, init_logger

logger = logging.getLogger(__name__)

class HelpLogger(EnvExperiment):
    """help logger

    Why doesn't this produce a message in the artiq_gui Log?
    """

    def build(self):

    def prepare(self):
        logger.info("prepare()")

    def run(self):
        logger.info("run()")

    def analyze(self):
        logger.info("analyze()")

if __name__ == "__main__":
    from artiq.frontend.artiq_run import run
    run()