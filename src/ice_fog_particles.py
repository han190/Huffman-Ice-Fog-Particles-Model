import numpy as np
import matplotlib.pyplot as plt


class iceFogParticles:

    Tm = 273.15
    D0 = 0.106144
    R = 8.31446261815324e7
    M = 18.01528
    K = 2.09e3

    def __init__(self, v0, a, b, Ti, T0, f):
        self.v0 = v0
        self.a = a
        self.b = b
        self.Ti = Ti
        self.T0 = T0
        self.f = f

    def T(self, t):
        c = self.v0 / self.b
        deltaT = Ti - T0
        bottom = self.a * self.b * deltaT * np.log(c * t + 1.) + 1.
        return T0 + deltaT / bottom

    def dTdt(self, t):
        c = self.v0 / self.b
        deltaT = Ti - T0
        bottom = self.a * self.b * deltaT * np.log(c * t + 1.) + 1.
        top = self.a * self.b * self.c * deltaT**2
        return - top / ((c * t + 1.) * bottom**2)

    def Ew(self, t):
        c = [8.72679635, 27.84317001, -17.82965432]
        temp = T(self, t)
        termOne = c[1] * (1. - Tm / temp)
        termTwo = c[2] * np.log10(temp/Tm)
        return np.exp(c[0] + termOne + termTwo)

    def Ei(self, t):
        c = [8.9498907, -9.04755449, -36.00978668]
        temp = T(self, t)
        termOne = c[1] * (Tm / temp - 1.)
        termTwo = c[2] * np.log10(Tm/temp)
        return np.exp(c[0] + termOne + termTwo)

    def Es(self, t):
        temp = T(self, t)
        ew = Ew(self, t)
        ei = Ei(self, t)

        if temp > Tm:
            return ew
        else:
            return e_w + (T_m - temp) / 40. * (e_i - e_w)



