"""Library of structural statics-related functions"""
__author__ = "Ben Fisher"
__version__ = "0.0.22"
__license__ = "GPL"
__credits__ = ["Ben Fisher"]
__status__ = "Development"

import numpy as np

_STATIONS = 100

class Unit_Beam:
    def __init__(self):
        unitL = np.arange(0,1+1/_STATIONS,1/_STATIONS) 
        self.unitL = unitL
        self.R1 = 0                             # Reaction at left support
        self.R2 = 0                             # Reaction at right support
        self.W = np.zeros((1,_STATIONS))        # Distributed load array
        self.V = np.array([])                   # Shear vector
        self.M = np.array([])                   # Moment vector
        self.Rot = np.array([])                 # Rotation vector
        self.Def = np.array([])                 # Deflection vector

class Beam(Unit_Beam):
    def __init__(self, length):
        super().__init__()
        self.length = length
        self.x = length * self.unitL
        self.stations = self.x.size

class Load_DistrZero():
    def __init__(self, mag = 0):
        self.W = np.full([1,_STATIONS],mag)