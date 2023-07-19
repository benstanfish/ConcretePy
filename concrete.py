import materials
from math import sqrt, copysign, abs

def get_a(bw, Astl, concrete: materials.ConcreteMaterial, steel: materials.RebarMaterial):
    try:
        return steel.fy * Astl / (0.85 * concrete.fc * bw)
    except ZeroDivisionError:
        return 0

def get_c_from_a(a, concrete: materials.ConcreteMaterial):
    try:
        return a / concrete.b1
    except ZeroDivisionError:
        return 0

def get_c_from_Z(z, d, concrete: materials.ConcreteMaterial, steel: materials.RebarMaterial):
    try:
        return d / (1 - z * steel.ey / concrete.ecu)
    except ZeroDivisionError:
        return 0

def get_c_from_d(strain:float , d:float , concrete: materials.ConcreteMaterial):
    try:
        return d / (1 - strain / concrete.ecu)
    except ZeroDivisionError:
        return 0

def get_es(c, d, concrete: materials.ConcreteMaterial):
    try:
        return concrete.ecu*(1-d/c)
    except ZeroDivisionError:
        return 0

def get_fs(strain: float, steel: materials.RebarMaterial):
    if abs(strain/steel.ey) < 1:
        return copysign(strain*steel.Es,strain)
    else:
        return copysign(steel.fy,strain)

def geometric_sequence(n, initial, common_ratio):
    return initial * common_ratio ^ (n - 1)
