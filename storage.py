



def get_Vc(bw, d, concrete: materials.ConcreteMaterial):
    return 2 * concrete.lam * sqrt(concrete.fc) * bw * d

def get_Vc_with_axial(bw, d, Nu, Ag, concrete: materials.ConcreteMaterial, isSignificantTension: bool = False):
    # Eq. (22.5.6.1) or (22.5.7.1)
    denominator = 2000
    if isSignificantTension == True:
        denominator = 500
    return max((1 + Nu/(denominator * Ag)) * get_Vc(bw, d, concrete),0)

def get_Av_s(Vu, Vc, d, steel: materials.SteelMaterial, isSeisShear: bool = False):
    try:
        phiV = 0.75
        if isSeisShear == True:
            phiV = 0.6
        return (Vu - phiV * Vc) / (phiV * steel.fy * d)
    except ZeroDivisionError:
        return 0

def get_ac(hw: float, lw: float):
    # Alpha_c per ACI 318 Sec. 18.10.4.1
    try:
        if hw/lw <= 1.5:
            return 3
        elif hw/lw >= 2.0:
            return 2
        else:
            return 3 - 2 * (hw/lw - 1.5)
    except ZeroDivisionError:
        return 0