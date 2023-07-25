"""Library of functions to create PM coordinates for creating concrete PM diagrams"""
__version__ = "0.0.4"

print('concrete.py <version {}> successfully imported'.format(__version__))

import matplotlib.pyplot as plt
import numpy as np

import materials as mat

comp_c = 1e6  # Pure compression condition: c = infinity
tens_c = 0  # Pure tension condition: c = 0

def geometric_sequence(n, initial, common_ratio):
    return initial * common_ratio ^ (n - 1)

def equal_layer_distances(layer_count, 
                          bar_diameter, 
                          clear_cover, 
                          total_member_depth):
    return np.linspace(clear_cover + bar_diameter/2, 
                       total_member_depth - clear_cover - bar_diameter/2, 
                       layer_count)

def reverse_layers(total_member_depth, 
                   layer_distances):
    layers = layer_distances.copy()
    return np.flip(total_member_depth - layers)

def max_axial(gross_area, layer_areas, 
             concrete: mat.ConcreteMaterial, 
             rebar: mat.RebarMaterial, 
             isTensionCase: bool = False):
    if isTensionCase == False:
        return (0.85 * concrete.fc) * (gross_area - np.sum(layer_areas)) + rebar.fy * np.sum(layer_areas)
    else:
        return np.sum(layer_areas) * rebar.fy * -1

#================================================================================
#    Using similar triangles, one can relate the concrete strain (ecu) with that
#    in the rebar at a depth "d" and the neutral axis depth "c". Using the
#    coefficient Z as the ratio of steel strain at "d" to the yield strain "ey"
#================================================================================

def max_z(concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return concrete.ecu / (rebar.fy / rebar.Es)  # Occurs at maximum (pure) compression, c = inf

def min_z(rebar: mat.RebarMaterial):
    return rebar.eu / (rebar.fy / rebar.Es)  # Occurs at maximum (pure) tension; c = 0

def comp_z(concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    # This is just a wrapper function
    return max_z(concrete, rebar)  

def tens_z(rebar: mat.RebarMaterial):
    # This is just a wrapper function
    return min_z(rebar)

def z_from_c(c, layer_distance, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if c == 0:
        return max_z(concrete, rebar)
    else:
        return (c - layer_distance) / c * (concrete.ecu / rebar.ey)

def c_from_z(z, layer_distance, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if z * rebar.ey == concrete.ecu:
        return comp_c
    else:
        return max(layer_distance / (1 - (z * rebar.ey / concrete.ecu)), 0)

def c_from_strain(layer_strain, layer_distance, concrete: mat.ConcreteMaterial):
    if layer_strain == concrete.ecu:
        return comp_c
    else:
        return layer_distance / (1 - (layer_strain / concrete.ecu))

#================================================================================
#    
#                Formulae for calculating individual P-M points
#    
#================================================================================    

def layer_strain(layer_distance, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if c == 0:
        return rebar.eu
    else:
        return (c - layer_distance) * concrete.ecu / c

def layer_stress(layer_distance, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    strain = layer_strain(layer_distance, c, concrete, rebar)
    if abs(strain) >= rebar.ey:
        return rebar.fy * np.sign(strain)
    else:
        return strain * rebar.Es
    
def layer_force(layer_area, layer_distance, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if layer_distance >= c:
        return layer_stress(layer_distance, c, concrete, rebar) * layer_area
    else:
        return (layer_stress(layer_distance, c, concrete, rebar) - 0.85 * concrete.fc) * layer_area
    
def sum_layer_forces(layer_areas, layer_distances, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    # This function only considers the steel contributions.
    force = 0
    for i in range(layer_distances.shape[0]):
        force += layer_force(layer_areas[i], layer_distances[i], c, concrete, rebar)
    return force

def sum_layer_moments(layer_areas, layer_distances, c, h, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    # This function only considers the steel contributions.
    moment = 0
    for i in range(layer_distances.shape[0]):
        moment += layer_force(layer_areas[i], layer_distances[i], c, concrete, rebar) * (h/2 - layer_distances[i])
    return moment

def sum_forces(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    steel_force = sum_layer_forces(layer_areas, layer_distances, c, concrete, rebar)
    concrete_force = 0.85 * bw * (min(c * concrete.beta1, h)) * concrete.fc
    return steel_force + concrete_force

def sum_moments(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    steel_moment = sum_layer_moments(layer_areas, layer_distances, c, h, concrete, rebar)
    concrete_force = 0.85 * bw * min(c * concrete.beta1, h) * concrete.fc
    concrete_moment = concrete_force * (h - (min(c * concrete.beta1, h))) / 2
    return steel_moment + concrete_moment       
    
def pm_points(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    P = sum_forces(c, bw, h, layer_distances, layer_areas, concrete, rebar)
    M = sum_moments(c, bw, h, layer_distances, layer_areas, concrete, rebar)
    return P, M

#================================================================================
#    
#             Additional Formulae for Calculating Specific P-M Points
#    
#================================================================================   

def z_at_p(p_goal, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    
    min_tolerance = 1e-12
    max_iterations = 50
    keep_running = True
    n = 0

    z_min = tens_z(rebar)
    z_max = comp_z(concrete, rebar)
    d = max(layer_distances)
    tolerance = (z_max - z_min) / 2
    
    while keep_running == True:
        c_min = c_from_z(z_min, d, concrete, rebar)
        c_max = c_from_z(z_max, d, concrete, rebar)
        p_min = sum_forces(c_min, bw, h, layer_distances, layer_areas, concrete, rebar) - p_goal
        p_max = sum_forces(c_max, bw, h, layer_distances, layer_areas, concrete, rebar) - p_goal
        
        z = (z_min + z_max) / 2
        c = c_from_z(z, d, concrete, rebar)
        p = sum_forces(c, bw, h, layer_distances, layer_areas, concrete, rebar) - p_goal
        
        if np.sign(p) == 0:
            return z
        elif np.sign(p) == np.sign(p_min):
            z_min = z
        else:
            z_max = z
        tolerance = (z_max - z_min) / 2
        n += 1
        if (n == max_iterations) | (tolerance <= min_tolerance):
            keep_running = False
    return z

def z_at_pure_m(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return z_at_p(0, bw, h, layer_distances, layer_areas, concrete, rebar)









def createCList(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    """Create list of 'c' values to be used for points on the PM curve."""
    #================================================================================
    #   Creates a list of "c" values for two sides of the PM diagram for columns
    #   and walls. By default it generates 25 points:
    #   1 and 25 - maximum compression
    #   2 and 24 - half way between Pmax and Z = 0
    #   3 and 23 - Z = 0 (es = 0)
    #   4 and 22 - Po, axial compression at 80% of Pmax
    #   5 and 21 - Z = -0.5 (es = 50% fy)
    #   6 and 20 - balanced failure (Z = -1)
    #   7, 19 - Half way between comp and tension control limits
    #   8 and 18 - es = -0.005 (tension control limit)
    #   9-11, 15-17 - three points from tens control to pure moment separated by the
    #                 geometric sequence (to give a good distribution)
    #   12, 14 - pure moment
    #   15 - pure tension
    #
    #   The list is a numpy array, that includes the full 360° rotation. The first
    #   half of the list is negative flexure, the back half is positive flexure.
    #   These two halves are generated as two separated arrays, sorted in opposite
    #   order, then combined into one. This ensures that the PM coordinates are 
    #   generated in the correct order (even if its different from the listing above)
    #================================================================================
    
    # First half of diagram
    
    c_0 = cAtMaxComp()  # Point 1
    c_3 = c_from_z(0, max(layer_distances), concrete, rebar)  # Point 3
    
    Pmax = max_axial(bw * h, layer_areas, concrete, rebar, isTensionCase=False)
    Zmax = max_z()
    Zmin = min_z()
    z_Po = zFromP(Zmax, Zmin, 0.8 * Pmax, bw, h, layer_distances, layer_areas, concrete, rebar)
    c_4 = c_from_z(z_Po, max(layer_distances), concrete, rebar)
    c_5 = c_from_z(-0.5, max(layer_distances), concrete, rebar) 
    c_6 = c_from_z(-1, max(layer_distances), concrete, rebar)
    c_8 = cFromStrain(-0.005, max(layer_distances), concrete, rebar)
    c_7 = (c_6 + c_8)/2
    Zm = zAtPureM(bw, h, layer_distances, layer_areas, concrete, rebar)
    c_12 = c_from_z(Zm, max(layer_distances), concrete, rebar)
    
    
    
    c_15 = tens_s()
    