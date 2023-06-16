"""Library of concrete-related functions based on the ACI 318"""
# Created by Ben Fisher, 2023
# Version 1.0

from math import sqrt, copysign
import csv

DIAMS = {
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

AREAS = {
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

WEIGHTS = {
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

DEFAULT_ES = 29000000.0    # Steel modulus of elasticity (Es) in psi per Sec. 20.2.2.2
DEFAULT_FY = 60000.0       # Default steel yeild strength in psi (Grade 60 bar)
DEFAULT_ECU = 0.003     # Default limit concrete strain. Negative strain == tension.

""" Arguments used:

    aDist (float): the "a" distance (in)
    Ag (float): gross cross-sectional area (in^2)
    Astl (float): area of steel at a given layer (in^2)
    Av (float): area of transverse steel (in^2)
    betaOne (float): value of beta1 per Table 22.2.2.4.3 (unitless)
    bw (float): member width (in)
    cDist (float): the "c" neutral axis distance (in)
    coords: list of [x, y] lists of coordinates of the polygon's vertices
    dDist (float): the "d" distance from compression fiber to steel layer (in)
    ecu (float): limit concrete compression strain, typically 0.003 (in/in)
    Es (float): modulus of elasticity of steel (psi)
    fc (float): 28-day concrete compression strength (psi)
    fy (float): steel yield strength (psi)
    fyt (float): transverse steel yield strength (psi). Defaults to DEFAULT_FY.
    isEShear (bool): is false unless shear from case described in Sec. 21.2.4.1
    isSignificantTension (bool): "judgement call" as to whether tension is "significant"
    lam (float): lightweight concrete factor, lambda (unitless)
    Nu (float): axial force (lbf), where compression is positive, tension is negative
    spacing (float): spacing of transverse steel (in)
    strain (float): steel strain (in/in) at layer "d"
    Vc (float): concrete shear strength (lbf)
    Vs (float): nominal transverse steel shear strenght (lbf) per Eq. (23.5.10.5.3)
    Vu (float): nominal shear demand (lbf)
    wc (float): concrete unit weight (pcf)
"""

def beta1(fc):
    """Calculate the value of beta1 per Table 22.2.2.4.3
    
    fc: 28-day concrete compression strength in psi"""
    # If ksi units are provided, this function converts to psi
    if fc < 10:
        fc *= 1000
    if fc <= 4000:
        return 0.85
    elif fc >= 8000:
        return 0.65
    else:
        return 0.85 - 0.05*(fc - 4000)/1000

def cDist(aDist, betaOne):
    """Returns neutral axis distance c per Eq. (22.2.2.4.1)
    
    aDist: the "a" distance (in)
    betaOne: value of beta1 per Table 22.2.2.4.3 (unitless)"""
    return aDist/betaOne

def aDist(fc, bw, Astl, fy = DEFAULT_FY):
    """Calculates the "a" distance of beam
    
    fc: 28-day concrete compression strength (psi)
    bw: member width (in)
    Astl: area of steel at a given layer (in^2)
    fy: steel yield strength (psi)"""
    # The factor 0.85f'c is the limit concrete stress per Sec. 22.2.2.4.1.
    return  fy * Astl / (0.85 * fc * bw)

def Ec(fc, wc: float = None):
    """Returns concrete elastic modulus (Ec) per Eq. (19.2.2.1.b), unless wc provided, then Eq. (19.2.2.1.a)
    
    fc: 28-day concrete compression strength (psi)
    wc: concrete unit weight (pcf), should be between 90 and 160 pcf
    """
    # Note: wc = 143.959593 pcf makes both equations approximately equal.
    # If ksi units are provided, this function converts to psi
    if fc < 10:
        fc *= 1000
    if wc is None:
        return 57000*sqrt(fc)
    else:
        return (wc**1.5)*33*sqrt(fc)

def ruptureModulus(fc, lam: float = 1):
    """Returns modulus of rupture per Eq. (19.2.3.1). Lambda per Table 19.2.4.2, defaults to lambda = 1 (NWC).

    fc: 28-day concrete compression strength (psi)
    lam: lightweight concrete factor, lambda (unitless)"""
    # If ksi units are provided, this function converts to psi
    if fc < 10:
        fc *= 1000
    return 7.5*lam*sqrt(fc)

def steelStrain(dDist, cDist, ecu: float = DEFAULT_ECU):
    """Returns the strain at layer "d" assuming similar triangles. Compression treated as positive.
    
    cDist: the "c" neutral axis distance (in)
    dDist: the "d" distance from compression fiber to steel layer (in)
    ecu = limit concrete strain of 0.003, per Sec. 22.2.2.1"""
    return ecu * (1 - dDist/cDist)

def yieldStrain(fy: float = DEFAULT_FY, Es: float = DEFAULT_ES):
    """Returns the yield strain (ey) for a given steel strength (psi) and elastic modulus (psi)

    fy: steel yield strength (psi)
    Es: modulus of elasticity of steel (psi)"""
    return fy/Es

def steelStress(strain, Es: float = DEFAULT_ES, fy: float = DEFAULT_FY):
    """Returns steel stress (fs) based on the steel strain. Refer to Sec. R20.2.2.1
    
    strain: steel strain (in/in) at layer "d"
    Es: modulus of elasticity of steel (psi)
    fy: steel yield strength (psi)
    """
    # Note that this case applies -ey < fyain < ey. Therefore steel in compression
    # must have a negative sign to get the correct stress sign
    ey = fy/Es
    if abs(strain/ey) < 1:
        return copysign(strain*Es,strain)
    else:
        return copysign(fy,strain)

def getVc(fc: float, bw: float, dDist: float, lam: float = 1):
    """Returns the shear strength Vc (lbf) for nonprestressed 
    members w/o axial forcce per Eq. (22.5.5.1)

    fc: 28-day concrete compression strength (psi)
    bw: member width (in)
    dDist: the "d" distance from compression fiber to steel layer (in)"""
    return 2*lam*sqrt(fc)*bw*dDist

def getVc_WithAxial(
    fc: float, bw: float, dDist: float, 
    Nu: float, Ag: float, lam: float = 1, isSignificantTension:bool = False):
    """Returns the shear strength Vc (lbf) for nonprestressed members WITH axial force per Eq. (22.5.6.1) or (22.5.7.1) if 
    
    Args:
        fc: 28-day concrete compression strength (psi)
        bw: member width (in)
        dDist: the "d" distance from compression fiber to steel layer (in)
        Nu: axial force (lbf), where compression is positive, tension is negative
        Ag: gross cross-sectional area (in^2)"""
    denom = 2000
    if isSignificantTension == True:
        denom = 500
    Vc = 2*(1+Nu/(denom*Ag))*lam*sqrt(fc)*bw*dDist
    return max(Vc,0)

def getVs(Av, spacing, dDist, fyt: float = DEFAULT_FY):
    """Returns the nominal shear strength per Eq. (23.5.10.5.3)

    Args:
        Av (float): area of transverse steel (in^2)
        spacing (float): spacing of transverse steel (in)
        dDist (float): the "d" distance from compression fiber to steel layer (in)
        fyt (float, optional): transverse steel yield strength (psi). Defaults to DEFAULT_FY.
    """
    return Av*fyt*dDist/spacing

def getAv_s(Vu, Vc, dDist, fyt: float = DEFAULT_FY, isEShear: bool = False):
    """Returns Av/s (in^2/in) per Eq. (R22.5.10.5) for direct shear

    Args:
        Vu (float): nominal shear demand (lbf)
        Vc (float): concrete shear strength (lbf)
        dDist (float): the "d" distance from compression fiber to steel layer (in)
        fyt (float): transverse steel yield strength (psi). Defaults to DEFAULT_FY.
        isEShear (bool): is false unless shear from case described in Sec. 21.2.4.1
    """
    # Note that this function assumes phi = 0.75 per Table 21.2.1(b), unless isSeismicShear
    # is True, in which case phi = 0.60 per Sec. 21.2.4.1
    phiV = 0.75
    if isEShear == True:
        phiV = 0.6
    return (Vu-phiV*Vc)/(phiV*fyt*dDist)

def getPhiVn_1W(Vs, Vc, isEShear: bool = False):
    """Returns phi*Vn for one-way shear, based Eq. (22.5.1.1) and Ch. 21 phi factors.

    Args:
        Vs (float): nominal transverse steel shear strenght (lbf) per Eq. (23.5.10.5.3)
        Vc (float): nominal one-way concrete shear strength (lbf)
        isEShear (bool, optional): _description_. Defaults to False.
    """
    phiV = 0.75
    if isEShear == True:
        phiV = 0.6
    return phiV*(Vc+Vs)







def area(coords):
    """Return the area of a closed, simple ("non-self-intersecting") 
    polygon based on x, y pairs; uses the trapezoid formula.

    Args:
        coords: list of [x, y] lists of coordinates of the polygon's vertices
    """
    # https://en.wikipedia.org/wiki/Shoelace_formula#Trapezoid_formula
    # Example: coords = [[-2, -2], [11, 2], [9, 7], [4, 10]] --> 75.5
    arr = coords.copy()     # Deep copy the coordiante array for the next step
    arr.append(coords[0])   # Copy the first indice to the end for later operations
    A = 0
    for i in range(0,len(arr)-1):
        A += 0.5*(arr[i][0]*arr[i+1][1]-arr[i+1][0]*arr[i][1])
    return abs(A)

def centroid(coords):
    """Return the centroid of a closed, simple ("non-self-intersecting") 
    polygon based on x, y pairs; uses a variation of the trapezoid formula.

    Args:
        coords: list of [x, y] lists of coordinates of the polygon's vertices
    """
    # Example: coords = [[-2, -2], [11, 2], [9, 7], [4, 10]] --> 75.5
    # https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    arr = coords.copy()     # Deep copy the coordiante array for the next step
    arr.append(coords[0])   # Copy the first indice to the end for later operations
    cgX = 0
    cgY = 0
    area = 0
    for i in range(0,len(arr)-1):
        area += 0.5*(arr[i][0] * arr[i+1][1] - arr[i+1][0] * arr[i][1])
    area = abs(area)
    for i in range(0,len(arr)-1):
        cgX += (arr[i][0]+arr[i+1][0])*(arr[i][0]*arr[i+1][1]-arr[i+1][0]*arr[i][1])/6/area
        cgY += (arr[i][1]+arr[i+1][1])*(arr[i][0]*arr[i+1][1]-arr[i+1][0]*arr[i][1])/6/area
    return [cgX, cgY, area]


# Further Research: https://www.spatialanalysisonline.com/HTML/centroids_and_centers.htm




# with open('points.csv', newline='') as csvfile:
#     myReader = csv.reader(csvfile, delimiter=" ", quotechar="|")
#     for row in myReader:
#         print(', '.join(row))


# with open('points.csv', "r", newline='\n') as csvfile:
#     myReader = csv.reader(csvfile, delimiter=' ')
#     arr = []
#     for row in myReader:
#         arr.append(row)
#     print(arr)

def readlines():
    with open('points.csv', 'r') as data:
        reader = csv.reader(data)
        for row in reader:
            yield [ float(i) for i in row ]

# for i in read_lines():
#     print(i)

# to get a list, instead of a generator, use
xy = list(readlines())

arr = centroid(xy)
print(arr[2])