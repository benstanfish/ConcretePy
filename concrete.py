"""Library of functions to create PM coordinates for creating concrete PM diagrams"""

import matplotlib.pyplot as plt
import numpy as np

import materials as mat

def geometric_sequence(n, initial, common_ratio):
    return initial * common_ratio ^ (n - 1)

def layerDistances(n, db, cc, h):
    return np.linspace(cc + db/2, h-cc-db/2, n)

def maximumCompression(Ag, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return (0.85 * concrete.fc) * (Ag - np.sum(layer_areas)) + rebar.fy * np.sum(layer_areas)

def cFromZ(Z, d, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    try:
        return concrete.ecu/(concrete.ecu - Z * rebar.ey)*d
    except ZeroDivisionError:
        return 0
    
def cFromD(strain_at_d, d, concrete: mat.ConcreteMaterial):
    try:
        return d / (1 - strain_at_d / concrete.ecu)
    except ZeroDivisionError:
        return 0
    
def zFromC(c, d, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if c != 0:
        try:
            return concrete.ecu / rebar.ey * (1 - d / abs(c))
        except:
            return 0
    else:
        # c = 0 is error case; pass 250 which is well past the rupture strain
        return 250

def layerStrain(layer_distance, c, concrete: mat.ConcreteMaterial):
    try:
        return (c - layer_distance)/c * concrete.ecu
    except:
        return 0
    
def layerStress(layer_distance, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    strain = layerStrain(layer_distance, c, concrete)
    return min(rebar.fy * np.sign(strain), strain * rebar.Es)

def layerForce(layer_area, layer_distance, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if layer_distance > c:
        return layerStress(layer_distance, c, concrete, rebar) * layer_area
    else:
        return (layerStress(layer_distance, c, concrete, rebar) - 0.85 * concrete.fc) * layer_area
    
def sumLayerForces(layer_areas, layer_distances, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    force = 0
    for i in range(layer_areas.shape[0]):
        force += layerForce(layer_areas[i], layer_distances[i], c, concrete, rebar)
    return force

    
def sumLayerMoments(layer_areas, layer_distances, h, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    moment = 0
    for i in range(layer_areas.shape[0]):
        moment += layerForce(layer_areas[i], layer_distances[i], c, concrete, rebar) * (h/2 - layer_distances[i])
    return moment

def PMPoints(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    sum_Fs = 0
    sum_Ms = 0
    for i in range(layer_areas.shape[0]):
        sum_Fs += layerForce(layer_areas[i], layer_distances[i], c, concrete, rebar)
        sum_Ms += layerForce(layer_areas[i], layer_distances[i], c, concrete, rebar) * (h/2 - layer_distances[i])
    Cc = 0.85 * bw * (c * concrete.b1) * concrete.fc
    Mc = Cc * (h - (c * concrete.b1))/2
    P = Cc + sum_Fs
    M = Mc + sum_Ms
    return P, M

def ZfromP(Za, Zb, P, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    """Find Z (relating to ey) for a given P value. Solve using the Newton-Rapfson root finding method.
    Calculated Pi values are shifted down by P (such that P => y = 0) to solve for it as a root.
    The initial values of Za and Zb must result in P values above and below P argument, i.e., bound the desired Z value."""
    
    minTol = 0.000000000000001  # Value selected to push error just outside the data type smallest value
    maxIter = 500  # The maximum iterations does not usually control.
    keepRunning = True
    n = 0
    d = max(layer_distances)
    tol = (Zb - Za) / 2
    
    while keepRunning == True:
        ca = cFromZ(Za, d, concrete, rebar)
        Pa = PMPoints(ca, bw, h, layer_distances, layer_areas, concrete, rebar)[0] - P

        Zc = (Za + Zb) / 2
        print(str(n).rjust(2,"0"), Zc)  # This is a debug string
        cc = cFromZ(Zc, d, concrete, rebar)
        Pc = PMPoints(cc, bw, h, layer_distances, layer_areas, concrete, rebar)[0] - P 

        if abs(Pc) == 0:
            return Zc
        elif np.sign(Pc) ==  np.sign(Pa):
            Za = Zc
        else:
            Zb = Zc
        
        tol = (Zb - Za) / 2
        n += 1
        if n == maxIter or tol <= minTol:
            keepRunning = False
    return Zc

def ZAtPureM(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    """Calculate Z (relating to ey) where P = 0, i.e, the pure moment condition."""
    Za = 250  # Selected to be sufficiently high.
    Zb = 0
    P = 0
    return ZfromP(Za, Zb, P, bw, h, layer_distances, layer_areas, concrete, rebar)