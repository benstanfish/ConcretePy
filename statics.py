"""Library of structural statics-related functions"""
__author__ = "Ben Fisher"
__version__ = "0.0.32"
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
        self.W = np.zeros((1,1+_STATIONS))        # Distributed load array
        self.V = np.array([])                   # Shear vector
        self.M = np.array([])                   # Moment vector
        self.Rot = np.array([])                 # Rotation vector
        self.Def = np.array([])                 # Deflection vector

    def addDistr(self, distr_load):
        self.W = np.add(self.W, distr_load.W)

class Beam(Unit_Beam):
    def __init__(self, length):
        super().__init__()
        self.length = length
        self.x = length * self.unitL

class Distr_Load():
    """
    _summary_
    """
    def __init__(self, mag = 0, a = 0.0, b = 1.0, c = 0.0):
        """
        _summary_

        Args:
            mag (float, optional): Magnitude of load. Defaults to 0.0.
            a (float, optional): Distance to load from left support. Defaults to 0.0.
            b (float, optional): Length of loaded area. Defaults to 1.0, full length.
            c (float, optional): Distance from load to right support. Defaults to 0.0.
        """
        self.unitL = np.arange(0,1+1/_STATIONS,1/_STATIONS)
        self.W = np.copy(self.unitL)
        self.W = np.where((self.W >= a) & (self.W <= 1 - c),mag,0)
        
class SimpleUniform():
    def __init__(self, magnitude, length):
        unitL = np.arange(0,1+1/_STATIONS,1/_STATIONS)
        self.unitL = unitL
        self.x = unitL*length
        self.w = magnitude
        self.length = length
        self.R1 = self.R2 = magnitude * length / 2
        self.V = magnitude * (length/2 - self.x)