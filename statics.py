"""Library of structural statics-related functions"""
__author__ = "Ben Fisher"
__version__ = "0.0.36"
__license__ = "GPL"
__credits__ = ["Ben Fisher"]
__status__ = "Development"

import numpy as np

_STATIONS = 100

class Uniform():
    def __init__(self, w=1.0, length=1.0, E = 1.0, I = 1.0):
        unitL = np.arange(0,1+1/_STATIONS,1/_STATIONS)
        self.length = length
        self.unitL = unitL
        self.x = unitL * length
        self.w = w
        self.E = E
        self.I = I
        
        self.R1 = self.R2 = w * length / 2
        self.V = w * (length / 2 - self.x)
        self.M = w / 2 * self.x * (length - self.x)
        
        self.d = -( w / (24 * E * I) * self.x * (length**3 - 2 * length * self.x**2 + self.x**3))
        
        self.vmin = np.min(self.V)
        self.vmax = np.max(self.V)
        self.mmax = np.min(self.M)
        self.mmax = np.max(self.M)
        self.dmin = np.min(self.d)
        self.dmax = np.max(self.d)
        
        self.xvmin = self.x[np.where(self.V == np.min(self.V))[0][0]]
        self.xvmax = self.x[np.where(self.V == np.max(self.V))[0][0]]
        self.xmmin = self.x[np.where(self.M == np.min(self.M))[0][0]]
        self.xmmax = self.x[np.where(self.M == np.max(self.M))[0][0]]
        self.xdmin = self.x[np.where(self.d == np.min(self.d))[0][0]]
        self.xdmax = self.x[np.where(self.d == np.max(self.d))[0][0]]