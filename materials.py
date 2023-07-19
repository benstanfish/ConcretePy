from math import sqrt
from typing import List, Union, Any, Dict

# default properties for US materials using ACI 318 and imperial units

# reinforcement values per ACI 318, Appendix A
bar_numbers = [3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 18]
bar_diameters = {
    3: 0.375,
    4: 0.500,
    5: 0.625,
    6: 0.750,
    7: 0.875,
    8: 1.000,
    9: 1.128,
    10: 1.270,
    11: 1.410,
    14: 1.693,
    18: 2.257
}
bar_areas = {
    3: 0.11,
    4: 0.20,
    5: 0.31,
    6: 0.44,
    7: 0.60,
    8: 0.79,
    9: 1.00,
    10: 1.27,
    11: 1.56,
    14: 2.25,
    18: 4.00
}
bar_weights = {
    3: 0.376,
    4: 0.668,
    5: 1.043,
    6: 1.503,
    7: 2.044,
    8: 2.67,
    9: 3.4,
    10: 4.303,
    11: 5.313,
    14: 7.65,
    18: 13.6
}

class ConcreteMaterial:
    def __init__(self, fc, lam=1.0):
        self.fc = fc  # Min 28-day strength (psi)
        self.lam = lam  # LWC factor lambda
        self.ecu = 0.003  # Crushing strain
        self.Ec = 57000 * sqrt(self.fc)  # Elastic modulus (psi)
        self.fr = 7.5 * lam * sqrt(fc)  # Modulus of rupture (psi)
        self.b1 = self.beta1

    @property
    def beta1(self):
        if self.fc <= 4000:
            return 0.85
        elif self.fc >= 8000:
            return 0.65
        else:
            return 0.85 - 0.05 * (self.fc - 4000) / 1000


class RebarMaterial:
    def __init__(self, fy = 60000):
        self.fy = fy  # Mim steel yield stress (psi)
        self.Es = 29000000  # Elastic modulus (psi)
        self.ey = self.fy / self.Es  # Steel yield strain
