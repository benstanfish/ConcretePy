'''Library based on 2018 ACI 318, Imperial Units'''
# Version 1.0

from math import sqrt, copysign
from typing import overload

diams = {
    # Nominal linear diameter of rebar per Appendix A
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

areas = {
    # Nominal cross-sectional area of rebar per Appendix A
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

weights = {
    # Nominal linear weight of rebar per Appendix A
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

defaultEs = 29000000    # Steel modulus of elasticity (Es) in psi per Sec. 20.2.2.2
defaultFy = 60000       # Default steel yeild strength in psi (Grade 60 bar)
defaultEcu = -0.003     # Default limit concrete strain. Negative strain == compression.

def get_beta1(fprimec: float):
    '''Returns beta1 value per Table 22.2.2.4.3
    
    fprimec = 28-day concrete compression strength in psi'''
    # If ksi units are provided, this function converts to psi
    if fprimec < 10:
        fprimec *= 1000
    if fprimec <= 4000:
        return 0.85
    elif fprimec >= 8000:
        return 0.65
    else:
        return 0.85 - 0.05*(fprimec - 4000)/1000

def get_c(aDist: float, beta1: float):
    '''Returns neutral axis distance "c" per Eq. (22.2.2.4.1)
    
    aDist = "a" distance in inches
    beta1 = value per Table 22.2.2.4.3'''
    return aDist/beta1

def get_a(fprimec: float, bw: float, Astl: float, fy: float = defaultFy):
    '''Returns "a" distance, assuming fs >= fy
    
    fprimec = 28-day concrete compression strength in psi
    bw = width of member in inches
    steelArea = total area of steel
    steelStr = steel strenght in psi, defaults to 60,000 psi'''
    # The factor 0.85f'c is the limit concrete stress per Sec. 22.2.2.4.1.
    return Astl * fy / (0.85 * fprimec * bw)

def get_Ec(fprimec: float, wc: float = None):
    '''Returns concrete elastic modulus, Ec per Eq. (19.2.2.1.b) unless wc provided, then Eq. (19.2.2.1.a)
    
    fprimec = 28-day concrete compression strength in psi
    wc = concrete unit weight in lb/ft^3, between 90 and 160 pcf.
    
    Note that wc = 143.959593 pcf makes both equations approximately equal.
    '''
    # If ksi units are provided, this function converts to psi
    if fprimec < 10:
        fprimec *= 1000
    if wc is None:
        return 57000*sqrt(fprimec)
    else:
        return (wc**1.5)*33*sqrt(fprimec)

def get_fr(fprimec: float, lam: float = 1):
    '''Returns modulus of rupture per Eq. (19.2.3.1). Note lambda per Table 19.2.4.2, defaults to lambda = 1 (NWC).'''
        # If ksi units are provided, this function converts to psi
    if fprimec < 10:
        fprimec *= 1000
    return 7.5*lam*sqrt(fprimec)

def get_strain_at_d(dDist: float, cDist: float, ecu: float = defaultEcu):
    '''Returns the strain at layer "d" assuming similar triangles. Tension = positive strains.
    
    ecu = limit concrete strain of -0.003, per Sec. 22.2.2.1 (negative for compression).'''
    return ecu * (1 - dDist/cDist)

def get_fy(fy: float = defaultFy, Es:float = defaultEs):
    '''Returns the yield stress (fy) for a given steel strength (psi) and elastic modulus (psi)'''
    return fy/Es

def get_fs(steelStrain: float, Es: float = defaultEs, fy: float = defaultFy):
    '''Returns steel stress (fs) based on the steel strain. Refer to Sec. R20.2.2.1'''
    # Note that this case applies -ey < steelStrain < ey. Therefore steel in compression
    # must have a negative sign to get the correct stress sign
    ey = fy/Es
    if abs(steelStrain/ey) < 1:
        return copysign(steelStrain*Es,steelStrain)
    else:
        return copysign(fy,steelStrain)

def get_Vc(fprimec: float, bw: float, dDist: float, lam: float = 1):
    '''Returns the shear strength Vc (lbf) for nonprestressed 
    members w/o axial forcce per Eq. (22.5.5.1)

    fprimec = 28-day concrete compression strength (psi)
    bw = member width in inches
    dDist = distance to extreme tension steel layer in inches'''
    return 2*lam*sqrt(fprimec)*bw*dDist

def get_Vc_with_axial(
    fprimec: float, bw: float, dDist: float, 
    Nu: float, Ag: float, lam: float = 1, isSignificantTension:bool = False):
    '''Returns the shear strength Vc (lbf) for nonprestressed 
    members WITH axial force per Eq. (22.5.6.1) or (22.5.7.1) if 
    
    fprimec = 28-day concrete compression strength (psi)
    bw = member width in inches
    dDist = distance to extreme tension steel layer in inches
    Nu = axial force (lbf), where compression is positive, tension is negative
    Ag = gross cross-sectional area (in^2)'''
    denom = 2000
    if isSignificantTension == True:
        denom = 500
    Vc = 2*(1+Nu/(denom*Ag))*lam*sqrt(fprimec)*bw*dDist
    return max(Vc,0)

