"""Library of structural statics-related functions"""
__author__ = "Ben Fisher"
__version__ = "0.0.36"
__license__ = "GPL"
__credits__ = ["Ben Fisher"]
__status__ = "Development"

import numpy as np

station_count = 101

class Uniform():
    def __init__(self, w=1.0, length=1.0, E = 1.0, I = 1.0):
        unitL = np.arange(0,1+1/station_count,1/station_count)
        self.length = length
        self.unitL = unitL
        self.Xs = unitL * length
        self.w = w
        self.E = E
        self.I = I
        
        self.R1 = self.R2 = w * length / 2
        self.V = w * (length / 2 - self.Xs)
        self.M = w / 2 * self.Xs * (length - self.Xs)
        
        self.d = -( w / (24 * E * I) * self.Xs * (length**3 - 2 * length * self.Xs**2 + self.Xs**3))
        
        self.vmin = np.min(self.V)
        self.vmax = np.max(self.V)
        self.mmax = np.min(self.M)
        self.mmax = np.max(self.M)
        self.dmin = np.min(self.d)
        self.dmax = np.max(self.d)
        
        self.vmin_x = self.Xs[np.where(self.V == np.min(self.V))[0][0]]
        self.vmax_x = self.Xs[np.where(self.V == np.max(self.V))[0][0]]
        self.mmin_x = self.Xs[np.where(self.M == np.min(self.M))[0][0]]
        self.mmax_x = self.Xs[np.where(self.M == np.max(self.M))[0][0]]
        self.dmin_x = self.Xs[np.where(self.d == np.min(self.d))[0][0]]
        self.dmax_x = self.Xs[np.where(self.d == np.max(self.d))[0][0]]
        

class UniformDist_Pct():
    """_summary_
    """
    def __init__(self, mag = 1, length = 1, a = 0, c = 0):
        """_summary_

        Args:
            mag (int, optional): _description_. Defaults to 1.
            a (int, optional): _description_. Defaults to 0.
            c (int, optional): _description_. Defaults to 0.
        """
        self.b = length - a - c
        if self.b < 0:
            self.b = 0
        elif self.b > length:
            self.b = length
        
        self.Xs = np.linspace(0,1,station_count)
        
        

        