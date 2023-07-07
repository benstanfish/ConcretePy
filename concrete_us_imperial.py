import materials
from math import sqrt, copysign

def get_a(bw, Astl, concrete: materials.ConcreteMaterial,
          steel: materials.SteelMaterial):
    try:
        return steel.fy * Astl / (0.85 * concrete.fc * bw)
    except ZeroDivisionError:
        return 0

def get_c(a, concrete: materials.ConcreteMaterial):
    try:
        return a / concrete.b1
    except ZeroDivisionError:
        return 0

def get_es(c,d,concrete: materials.ConcreteMaterial):
    try:
        return concrete.ecu*(1-d/c)
    except ZeroDivisionError:
        return 0

def get_fs(strain: float, steel: materials.SteelMaterial):

