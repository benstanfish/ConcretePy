from math import sqrt

# default properties for US materials using ACI 318 and imperial units

# reinforcement values per ACI 318, Appendix A
bar_numbers = [3,4,5,6,7,8,9,10,11,14,18]
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

class Concrete_Material:
    def __int__(self, fc, lam=1.0):
        self.fc = fc                    # Min 28-day strength (psi)
        self.Ec = 57000*sqrt(self.fc)   # Elastic modulus (psi)
        self.fr = 7.5*lam*sqrt(fc)      # Modulus of rupture (psi)
        self.lam = lam                  # Lightweight concrete factor lambda
        self.ecu = 0.003                # Crushing strain
        
class Steel_Material:
    def __init__(self, fy):
        self.fy = fy                    # Mim steel yield stress (psi)
        self.Es = 29000000              # Elastic modulus (psi)
        self.ey = self.fy/self.Es       # Steel yield strain
        
    