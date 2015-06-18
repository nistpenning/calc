# vi: ts=4 sw=4

class TrapConfiguration():
    def __init__(self):
        self.Bz = 0
        self.kx = 0
        self.ky = 0
        self.kz = 0
        self.omega = 0
        self.theta = 0

    def changeDelta(self,delta):
        self.kx = -(0.5 + delta) * self.kz
        self.ky = -(0.5 - delta) * self.kz
