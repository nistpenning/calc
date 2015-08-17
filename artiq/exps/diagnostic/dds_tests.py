from artiq import *
from artiq.coredevice.dds import PHASE_MODE_CONTINUOUS
from artiq.coredevice.dds import PHASE_MODE_ABSOLUTE
from artiq.coredevice.dds import PHASE_MODE_TRACKING
import sys


class TestDDS(EnvExperiment):
    @kernel
    def mdelay(self, d):
        self.ttl1.pulse(d*0.25)
        delay(d*0.75)

    """DDS Test PHASE_MODE_ABSOLUTE"""
    def build(self):
        self.attr_device("core")
        self.attr_device("dds_bus")
        # to switch up which DDSs are used edit ddb.pyon
        self.attr_device("ddsA")
        self.attr_device("ddsB")
        self.attr_device("ddsC")
        self.attr_device("ttl0")
        self.attr_device("ttl1")

    @kernel
    def test5(self):
        """RF phase alignment with logic output"""
        # this works

        self.ddsA.set_phase_mode(PHASE_MODE_ABSOLUTE)
        while True:
            # mark start of experiment
            self.ttl0.pulse(2*us)
            self.ddsA.set(1*MHz, phase=0)
            delay(2.443*ms)  

    @kernel
    def test0(self):
        """absolute phase set at outset; do RF phase increment"""
        # This fails to change phase.

        self.ddsA.set_phase_mode(PHASE_MODE_ABSOLUTE)
        while True:
            # mark start of experiment
            self.ttl0.pulse(2*us)
            self.ddsA.set(2*MHz, phase=0)
            delay(2*us)
            self.ddsA.set(2*MHz, phase=0.5)
            delay(2*us)
            self.ddsA.set(2*MHz, phase=0)
            delay(2*us)
            self.ddsA.set(10*MHz, phase=0)
            delay(1*ms)
    @kernel
    def test4(self):
        """using ddsN.set(freq, phase_mode=PHASE_MODE_ABSOLUTE); do RF phase increment"""
        # This fails to change phase.

        while True:
            # mark start of experiment
            self.ttl0.pulse(2*us)
            self.ddsA.set(2*MHz, phase=0, phase_mode=PHASE_MODE_ABSOLUTE)
            self.mdelay(2*us)
            self.ddsA.set(2*MHz, phase=0.5, phase_mode=PHASE_MODE_ABSOLUTE)
            self.mdelay(2*us)
            self.ddsA.set(2*MHz, phase=0, phase_mode=PHASE_MODE_ABSOLUTE)
            self.mdelay(2*us)
            self.ddsA.set(10*MHz, phase=0, phase_mode=PHASE_MODE_ABSOLUTE)
            delay(1*ms)

    @kernel 
    def test1(self):
        # This fails to change phase.

        self.ddsA.set_phase_mode(PHASE_MODE_ABSOLUTE)
        self.ddsB.set_phase_mode(PHASE_MODE_ABSOLUTE)
        while True:
            self.ttl0.pulse(2*us)
            with self.dds_bus.batch:
                self.ddsA.set(1*MHz, phase=0)
                self.ddsB.set(1*MHz, phase=0)
            self.mdelay(2*us)
            with self.dds_bus.batch:
                self.ddsA.set(1*MHz, phase=0)
                self.ddsB.set(1*MHz, phase=0.5)
            self.mdelay(2*us)
            with self.dds_bus.batch:
                self.ddsA.set(1*MHz, phase=0)
                self.ddsB.set(1*MHz, phase=0)
            delay(1*ms)

    @kernel
    def test2(self):
        # PHASE_MODE_TRACKING works for ddsB comparing its phase before and after
        # frequency jump with reference ddsA. This fails however due to a slow 
        # phase drift between ddsA and ddsB. If line QQ is changed to 
        # phase_mode=PHASE_MODE_ABSOLUTE the drift goes away.

        while True:
            # force DDSs to be phase synchronized at outset of experiment
            with self.dds_bus.batch:
                self.ddsA.set(3*MHz, phase_mode=PHASE_MODE_ABSOLUTE, phase=0)
                self.ddsB.set(3*MHz, phase_mode=PHASE_MODE_ABSOLUTE, phase=0)

            # experiment start pulse on TTL0
            self.ttl0.pulse(2*us)
            self.ddsB.set(3*MHz, phase=0, phase_mode=PHASE_MODE_TRACKING) # QQ
            self.mdelay(2*us)
            # self.ddsB.set(10*MHz, phase=0, phase_mode=PHASE_MODE_TRACKING)
            # self.mdelay(2*us)
            # self.ddsB.set(3*MHz, phase=0, phase_mode=PHASE_MODE_TRACKING)
            # self.mdelay(2*us)
            # let core_uP catch up
            delay(1*ms)

    def test22(self):
        # Same level of functionality as test2

        while True:
            # force DDSs to be phase synchronized at outset of experiment
            with self.dds_bus.batch:
                self.ddsA.set(3*MHz, phase_mode=PHASE_MODE_ABSOLUTE, phase=0)
                self.ddsB.set(3*MHz, phase_mode=PHASE_MODE_ABSOLUTE, phase=0)

            # experiment start pulse on TTL0
            self.ttl0.pulse(2*us)
            self.ddsB.set(3*MHz, phase=0, phase_mode=PHASE_MODE_TRACKING)
            self.mdelay(2*us)
            self.ddsB.set(10*MHz, phase=0, phase_mode=PHASE_MODE_TRACKING)
            self.mdelay(2*us)
            self.ddsB.set(3*MHz, phase=0, phase_mode=PHASE_MODE_TRACKING)
            self.mdelay(2*us)
            delay(1*ms)

    @kernel
    def test3(self):

        # skip this test for now...

        # In this test mimic a spin-echo type pulse sequence.
        # ddsA is a reference tone
        # ddsB sets frequency of RF that is sent to atomic system.
        # ttl1 is a RF switch that exposes atoms to ddsB.
        # Suppose that the switch controlled by ttl1 has imprefect isolation and
        #   that we want to obtain additional decoupling from the atomic system
        #   by detuning during the spin-echo arm times. 

        self.ddsA.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.ddsB.set_phase_mode(PHASE_MODE_CONTINUOUS)
        tpi = 2*us
        f0 = 1*MHz
        df = 4*MHz # detuning to achieve better isolation from atomic system
        arm_t = 3*us
        dds_set_t = 0.5*us
        while True:
            # force DDSs to be phase synchronized at outset of experiment
            # this is a temporary override of the default phase mode
            with self.dds_bus.batch:
                self.ddsA.set(f0, phase_mode=PHASE_MODE_ABSOLUTE, phase=0)
                self.ddsB.set(f0+df, phase_mode=PHASE_MODE_ABSOLUTE, phase=0)
            # emit experiment start pulse on TTL0
            self.ttl0.pulse(2*us)

            #first pi/2
            self.ddsB.set(f0, phase=0)
            delay(dds_set_t) # wait for DDS to be programmed
            self.ttl1.pulse(tpi/2) 
            #arm 1
            self.ddsB.set(f0+df) # detune during arms 
            delay(arm_t)
            #middle pi
            self.ddsB.set(f0, phase=0.5)
            delay(dds_set_t) # wait for DDS to be programmed
            self.ttl1.pulse(tpi) # middle pi-pulse
            #arm 2
            self.ddsB.set(f0+df) # detune during arms 
            delay(arm_t)
            #final pi/2
            self.ddsB.set(f0)
            delay(dds_set_t) # wait for DDS to be programmed
            self.ttl1.pulse(tpi/2) # final pi/2 pulse

            # delay to let Core_uP catch up
            delay(1*ms)

    def input_choice(self):
        s = """Test DDSs in Penning lab.  
                Choose one of the following tests. 
                  5) PHASE_MODE_ABSOLUTE is RF aligned with logic?
                  0) PHASE_MODE_ABSOLUTE at outset
                  4) PHASE_MODE_ABSOLUTE at each set()
                  1) PHASE_MODE_ABSOLUTE with self.dds_bus.batch 
                  2) PHASE_MODE_TRACKING most basic 
                  22) PHASE_MODE_TRACKING fancy 
                  3) PHASE_MODE_CONTINUOUS with self.dds_bus.batch"""
        print(s)
        c = int(input("Choice: "))
        choices = {0: self.test0, 
            1: self.test1, 2: self.test2, 3: self.test3,
            4: self.test4, 5: self.test5}
        try:
            choices[c]()
        except KeyError:
            print("KeyError: Choose from [0,3]")
            sys.exit()

    def run(self):
        self.input_choice()
