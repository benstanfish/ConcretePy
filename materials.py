"""Library of material-related functions, based on ACI 318 and US customary units (inches, pounds)"""

from math import sqrt

#TODO: Global bar_* variables have been copied to RebarMaterial(), determine if they still need to be global here?
# Rebar sizes, diameters, areas, weights per ACI 318 Appendix A
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

# Material classes
class ConcreteMaterial:
    def __init__(self, fc, lam=1.0):
        self.fc = fc  # Minimum specified 28-day compressive strength (psi), user-specified.
        self.lam = lam  # Lambda factor for light-weight concrete (LWC) per ACI 318 Table 19.2.4.2.
        self.ecu = 0.003  # Maximum concrete compression strain per ACI 318 Sec. 22.2.2.1.
        self.Ec = 57000 * sqrt(self.fc)  # Elastic modulus (psi) per ACI 318 Eq. (19.2.2.1.b).
        self.fr = 7.5 * lam * sqrt(fc)  # Modulus of rupture (psi) per ACI 318 Eq. (19.2.3.1).
        self.b1 = self.beta1  # beta1 factor relating the "a" and "c" distances per ACI 318 Sec. 22.2.2.4.1.

    @property
    def beta1(self):
        """Calculate the beta1 factor per ACI 318 Eq. (22.2.2.4.1) and Table 22.2.2.4.3"""
        if self.fc <= 4000:
            return 0.85  # ACI 318 Table 22.2.2.4.3 case (a)
        elif self.fc >= 8000:
            return 0.65  # ACI 318 Table 22.2.2.4.3 case (c)
        else:
            return 0.85 - 0.05 * (self.fc - 4000) / 1000  # ACI 318 Table 22.2.2.4.3 case (b)

class SteelMaterial:
    def __init__(self, fy = 60000):
        self.fy = fy  # Minimum specified yeild-strength of rebar (psi), user-specified.
        self.Es = 29000000  # Modulus of elasticity for nonprestressed bars and wires per ACI 318 Sec. 20.2.2.2
        self.ey = self.fy / self.Es  # Yeild strain of the rebar, refer to ACI 318 Sec. R20.2.2.1
        
class RebarMaterial(SteelMaterial):
    def __init__(self, fy=60000):
        super().__init__(fy)
        self.bar_numbers = [3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 18]
        self.bar_diameters = {
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
        self.bar_areas = {
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
        self.bar_weights = {
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
        self.eu = self.ult_strain  # Rebar ultimate strain, defaults to ASTM A615 values.
        
    @property
    def ult_strain(self):
        """Ultimate strains associated with ASTM A615 material."""
        if self.fy <= 40000:
            return 0.155
        elif self.fy <= 60000:
            return 0.12
        elif self.fy <= 75000:
            return 0.07
        else:
            return 0.05  # For more information refer to ASCE 41-17, Sec. 10.3.3.1. 
                         # For historic rebar lower values may be more appropriate (e.g. 0.02, etc.)
    
    @property
    def reset_eu(self):
        """Reset the ultimate strain value to the original initialized value."""
        self.eu = self.ult_strain