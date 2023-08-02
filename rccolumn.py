"""Library of functions to create PM coordinates for creating concrete PM diagrams"""
__version__ = "0.0.4"
__author__ = "Ben Fisher"

import math

import matplotlib.pyplot as plt
import numpy as np

import rcmaterials as mat

CMAX = math.inf  # Pure compression condition: c = infinity
CCOMP = math.inf  # Pure compression condition: c = infinity
CMIN = 0  # Pure tension condition: c = 0
CTENS = 0  # Pure tension condition: c = 0

def equal_layer_distances(layer_count, bar_diameter, clear_cover, total_member_depth):
    return np.linspace(clear_cover + bar_diameter/2, total_member_depth - clear_cover - bar_diameter/2, layer_count)

def reverse_layers(total_member_depth, layer_distances):
    layers = layer_distances.copy()
    return np.flip(total_member_depth - layers)

def get_d(layer_distances):
    """Return d for the extreme tension steel layer"""
    return max(layer_distances)

def get_d_prime(layer_distances):
    """Return d' for the extreme compression steel layer"""
    return min(layer_distances)

def get_layer_areas(layer_bar_sizes, layer_bar_counts, rebar: mat.RebarMaterial):
    layers = layer_bar_counts.shape[0]
    layer_areas = np.zeros(layers)
    for i in range(layers):
        layer_areas[i] = rebar.bar_areas[layer_bar_sizes[i]] * layer_bar_counts[i]
    return layer_areas

def get_total_steel_area(layer_areas):
    return sum(layer_areas)

def column_reinforcement_ratio(concrete_gross_area, layer_areas):
    return sum(layer_areas) / concrete_gross_area

def max_axial(gross_area, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial, isTensionCase: bool = False):
    if isTensionCase == False:
        return (0.85 * concrete.fc) * (gross_area - np.sum(layer_areas)) + rebar.fy * np.sum(layer_areas)  # Po per Eq. (22.4.2.2)
    else:
        return np.sum(layer_areas) * rebar.fy * -1

def Po(gross_area, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return max_axial(gross_area, layer_areas, concrete, rebar, False)

def Pntmax(gross_area, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return max_axial(gross_area, layer_areas, concrete, rebar, True)

def capped_compression(max_compression, has_spirals: bool = False, is_ch_10_composite: bool = False):
    # max axial per ACI 318 Table 22.4.2.1
    # This function is a simplified version of Pnmax
    coeff = 0.80
    if has_spirals == True | is_ch_10_composite == True:
       coeff = 0.85
    return max_compression * coeff

def Pnmax(gross_area, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial, has_spirals: bool = False, is_ch_10_composite: bool = False):
    # Maximum compression strength, Pnmax, per ACI 318 Table 22.4.2.1
    coeff = 0.80
    if has_spirals == True | is_ch_10_composite == True:
       coeff = 0.85
    return coeff * max_axial(gross_area, layer_areas, concrete, rebar, isTensionCase = False)
    
    
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
        return min_z(rebar)  # Tension case
    else:
        return (c - layer_distance) / c * (concrete.ecu / rebar.ey)

def c_from_z(z, layer_distance, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if z * rebar.ey == concrete.ecu:
        return CCOMP
    else:
        return max(layer_distance / (1 - (z * rebar.ey / concrete.ecu)), 0)

def c_from_strain(layer_strain, layer_distance, concrete: mat.ConcreteMaterial):
    if layer_strain == concrete.ecu:
        return CCOMP
    else:
        return layer_distance / (1 - (layer_strain / concrete.ecu))

def cs_from_zs(zs, layer_distance, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    count = zs.shape[0]
    cs = np.zeros(count)
    for i in range(count):
        cs[i] = c_from_z(zs[i], layer_distance, concrete, rebar)
    return cs

def z_from_strain(strain, rebar: mat.RebarMaterial):
    return strain / rebar.ey

def zs_from_strains(strains, rebar: mat.RebarMaterial):
    """Batch create array of zs from array of strains"""
    count = strains.shape[0]
    zs = np.zeros(count)
    for i in range(count):
        zs[i] = strains[i] / rebar.ey
    return zs

def strain_from_c(c, layer_distance,  concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if c == 0:
        return rebar.eu
    else:
        return concrete.ecu * (1 - layer_distance / c)

def strain_from_z(z, rebar: mat.RebarMaterial):
    return z * rebar.ey

#================================================================================
#    
#                Formulae for calculating individual P-M points
#    
#================================================================================    

def layer_strain(layer_distance, c, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    if c == 0:
        return rebar.eu
    elif c == math.inf:
        # This case is necessary to prevent returning an "nan" error.
        return concrete.ecu
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

def sum_total_forces(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    steel_force = sum_layer_forces(layer_areas, layer_distances, c, concrete, rebar)
    concrete_force = 0.85 * bw * (min(c * concrete.beta1, h)) * concrete.fc
    return steel_force + concrete_force

def sum_total_moments(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    steel_moment = sum_layer_moments(layer_areas, layer_distances, c, h, concrete, rebar)
    concrete_force = 0.85 * bw * min(c * concrete.beta1, h) * concrete.fc
    concrete_moment = concrete_force * (h - (min(c * concrete.beta1, h))) / 2
    return steel_moment + concrete_moment       
    
def pm_points(c, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    P = sum_total_forces(c, bw, h, layer_distances, layer_areas, concrete, rebar)
    M = sum_total_moments(c, bw, h, layer_distances, layer_areas, concrete, rebar)
    d = max(layer_distances)
    strain_at_d = strain_from_c(c, d, concrete, rebar)
    return P, M, strain_at_d

def pm_from_cs(cs, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    count = cs.shape[0]
    P = np.zeros(count)
    M = np.zeros(count)
    strains = np.zeros(count)
    for i in range(count):
        P[i], M[i], strains[i] = pm_points(cs[i], bw, h, layer_distances, layer_areas, concrete, rebar)
    return P, M, strains

#================================================================================
#    
#             Additional Formulae for Calculating Specific P-M Points
#    
#================================================================================   

def get_z_at_p(p_goal, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    
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
        p_min = sum_total_forces(c_min, bw, h, layer_distances, layer_areas, concrete, rebar) - p_goal
        p_max = sum_total_forces(c_max, bw, h, layer_distances, layer_areas, concrete, rebar) - p_goal
        
        z = (z_min + z_max) / 2
        c = c_from_z(z, d, concrete, rebar)
        p = sum_total_forces(c, bw, h, layer_distances, layer_areas, concrete, rebar) - p_goal
        
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

def get_z_at_pure_m(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return get_z_at_p(0, bw, h, layer_distances, layer_areas, concrete, rebar)

def get_positive_zs(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, 
                       rebar: mat.RebarMaterial, 
                       points: int = 10, 
                       has_spirals: bool = False, 
                       is_ch_10_composite: bool = False):
        
    po = Po(bw * h, layer_areas, concrete, rebar)

    p_capped = Pnmax(bw * h, layer_areas, concrete, rebar, has_spirals, is_ch_10_composite)
    # z_at_p_capped = get_z_at_p(p_capped, bw, h, layer_distances, layer_areas, concrete, rebar)

    c_at_z_0 = c_from_z(0, max(layer_distances), concrete, rebar)
    p_at_z_0 = pm_points(c_at_z_0, bw, h, layer_distances, layer_areas, concrete, rebar)[0]
    
    a_range_step = (po - p_capped) / 5
    
    ps_a = np.arange(po, p_capped, -a_range_step)
    ps_b = np.linspace(p_capped, p_at_z_0, points - 4)    
    ps = np.hstack((ps_a, ps_b))
                           
    zs = np.zeros(ps.shape[0])
    for i in range(ps.shape[0]):
        zs[i] = get_z_at_p(ps[i], bw, h, layer_distances, layer_areas, concrete, rebar)
    return zs

def get_zs_from_zero_to_pure_m(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    """Create a list/array of z values from z = 0, through balanced failure, to tension control limit, to pure bending"""    
    cardinal_z_values = np.array([-0.125, -0.25, - 0.375, -0.5, -0.625, -0.75, -0.875, -1])
    cardinal_strains = np.arange(-0.0030, -0.006, -0.001)
    zs_from_cardinal_strains = zs_from_strains(cardinal_strains, rebar)
    
    z_Mmax = get_z_at_pure_m(bw, h, layer_distances, layer_areas, concrete, rebar)
    
    last_region_points = 5
    last_region_spaces = last_region_points - 1
    distance = zs_from_cardinal_strains[-1] - z_Mmax
    step_distance = distance / last_region_spaces
    
    last_region_zs = np.linspace(zs_from_cardinal_strains[-1] - step_distance, z_Mmax, last_region_points)
    return np.flip(np.sort(np.hstack((cardinal_z_values, zs_from_cardinal_strains, last_region_zs))))
    
def get_zs_in_tension_region(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    # TODO: need to figure out how to get Pnt in terms of z as z_from_p() doesn't appear to work in tension region.
    z_mmax = get_z_at_pure_m(bw, h, layer_distances, layer_areas, concrete, rebar)
    c_mmax = c_from_z(z_mmax, max(layer_distances), concrete, rebar)
    points = 9
    cs = np.linspace(0, c_mmax, points)
    zs = np.zeros(points)
    for i in range(points):
        zs[i] = z_from_c(cs[i], max(layer_distances), concrete, rebar)
    zs = np.insert(zs, 0, min_z(rebar) * 1000)  # This adds a really high point to ensure max Tension
    return zs

def get_half_zs(bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial, 
                points: int = 10, has_spirals: bool = False, is_ch_10_composite: bool = False):
    zs_positive = get_positive_zs(bw, h, layer_distances, layer_areas, concrete, 
                       rebar, points, has_spirals, is_ch_10_composite)
    zs_region_2 = get_zs_from_zero_to_pure_m(bw, h, layer_distances, layer_areas, concrete, rebar)
    zs_tension = get_zs_in_tension_region(bw, h, layer_distances, layer_areas, concrete, rebar)
    return np.sort(np.hstack((zs_positive, zs_region_2, zs_tension)))
    
def get_half_cs(zs, layer_distance, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    return cs_from_zs(zs, layer_distance, concrete, rebar)
    
def get_half_pm_points(cs, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial):
    P = np.zeros(cs.shape[0])
    M = np.zeros(cs.shape[0])
    strains = np.zeros(cs.shape[0])
    P, M, strains = pm_from_cs(cs, bw, h, layer_distances, layer_areas, concrete, rebar)
    return P, M, strains

def get_axial_moment_reduction_factor(strain_at_dt, rebar: mat.RebarMaterial, has_spirals: bool = False):
    """Return phi factor per ACI 318 Table 21.2.2 for axial, moment or axial + moment strength reduction factor"""
    phi_tens = 0.9
    phi_comp = 0.75 if has_spirals == True else 0.65
    coeff = 0.15 if has_spirals == True else 0.25
    return max(min(phi_comp + coeff * (strain_at_dt - -rebar.ey) / (-0.005 - -rebar.ey), phi_tens), phi_comp)

def get_multiple_axial_moment_reduction_factors(strains, rebar: mat.RebarMaterial, has_spirals: bool = False):
    count = strains.shape[0]
    phis = np.zeros(count)
    for i in range(count):
        phis[i] = get_axial_moment_reduction_factor(strains[i], rebar, has_spirals)
    return phis

def get_design_pm_points(cs, bw, h, layer_distances, layer_areas, concrete: mat.ConcreteMaterial, rebar: mat.RebarMaterial, has_spirals: bool = False, is_capped: bool = True):
    ps, ms, strains = get_half_pm_points(cs, bw, h, layer_distances, layer_areas, concrete, rebar)
    if is_capped == True:
        pnmax = Pnmax(bw * h, layer_areas, concrete, rebar)
        for i in range(ps.shape[0]):
            ps[i] = min(ps[i], pnmax)
    phis = get_multiple_axial_moment_reduction_factors(strains, rebar, has_spirals)
    design_ps = ps * phis
    design_ms = ms * phis
    return design_ps, design_ms













print(f'{__name__} <version {__version__}> successfully imported')