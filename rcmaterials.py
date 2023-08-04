"""Library of material-related functions, based on ACI 318 and US customary units (inches, pounds)"""
__version__ = "0.0.7"
__author__ = "Ben Fisher"


import numpy as np
import math

#TODO: Global bar_* variables have been copied to RebarMaterial(), determine if they still need to be global here?
# Rebar sizes, diameters, areas, weights per ACI 318 Appendix A
bar_numbers = [3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 18]

bar_diameters = {
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

bar_areas = {
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

bar_weights = {
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

# Material classes
class ConcreteMaterial:
    def __init__(self, fc, lam=1.0):
        self.fc = fc  # Minimum specified 28-day compressive strength (psi), user-specified.
        self.lam = lam  # Lambda factor for light-weight concrete (LWC) per ACI 318 Table 19.2.4.2.
        self.ecu = 0.003  # Maximum concrete compression strain per ACI 318 Sec. 22.2.2.1.
        self.Ec = 57000 * sqrt(self.fc)  # Elastic modulus (psi) per ACI 318 Eq. (19.2.2.1.b).
        self.fr = 7.5 * lam * sqrt(fc)  # Modulus of rupture (psi) per ACI 318 Eq. (19.2.3.1).
    
    @property
    def beta1(self):
        """Calculate the beta1 factor per ACI 318 Eq. (22.2.2.4.1) and Table 22.2.2.4.3"""
        if self.fc <= 4000:
            return 0.85  # ACI 318 Table 22.2.2.4.3 case (a)
        elif self.fc >= 8000:
            return 0.65  # ACI 318 Table 22.2.2.4.3 case (c)
        else:
            return 0.85 - 0.05 * (self.fc - 4000) / 1000  # ACI 318 Table 22.2.2.4.3 case (b)

class SteelMaterial:
    def __init__(self, fy = 60000):
        self.fy = fy  # Minimum specified yield-strength of rebar (psi), user-specified.
        self.Es = 29000000  # Modulus of elasticity for non-prestressed bars and wires per ACI 318 Sec. 20.2.2.2
        self.ey = self.fy / self.Es  # Yield strain of the rebar, refer to ACI 318 Sec. R20.2.2.1
        
class RebarMaterial(SteelMaterial):
    def __init__(self, fy=60000):
        super().__init__(fy)
        self.bar_numbers = [3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 18]
        self.bar_diameters = {
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
        self.bar_areas = {
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
        self.bar_weights = {
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
        self.eu = self.ult_strain  # Rebar ultimate strain, defaults to ASTM A615 values.
        
    @property
    def ult_strain(self):
        """Ultimate strains associated with ASTM A615 material."""
        if self.fy <= 40000:
            return -0.155
        elif self.fy <= 60000:
            return -0.12
        elif self.fy <= 75000:
            return -0.07
        else:
            return -0.05  # For more information refer to ASCE 41-17, Sec. 10.3.3.1. 
                         # For historic rebar lower values may be more appropriate (e.g. 0.02, etc.)
    
    @property
    def reset_eu(self):
        """Reset the ultimate strain value to the original initialized value."""
        self.eu = self.ult_strain
        
class RebarSet(RebarMaterial):
    def __init__(self, fy=60000):
        super().__init__(fy)
    
        
class ConcreteColumn(ConcreteMaterial):
    def __init__(self, widths_and_heights):
        super().__init__(fc, lam)
        self._cross_section = np.array(widths_and_heights)
        self._region_count = self.get_region_count(widths_and_heights)
        self.is_valid = self.check_if_valid(self.widths, self.heights)
    
    @property
    def cross_section(self):
        return self._cross_section

    @cross_section.setter
    def cross_section(self, new_widths_and_heights):
        self._cross_section = np.array(new_widths_and_heights)

    @property
    def widths(self):
        return self._cross_section[:,0]

    @property
    def heights(self):
        return self._cross_section[:,1]

    @property
    def gross_area(self):
        if self.is_valid:
            area = 0
            for i in range(self.heights.shape[0]):
                area += self.widths[i] * self.heights[i]
            return area
        else:
            return f"ERROR: cannot calculate area because there are {self.widths.shape[0]} width dimensions, and {self.heights.shape[0]} height dimensions."

    @property
    def gross_centroid(self):
        if self.is_valid:
            total_area = self.gross_area
            total_moment = 0
            for i in range(self.heights.shape[0]):
                total_moment += (self.widths[i] * self.heights[i] * (self.sum_cumulative(self.heights, False)[i] + self.heights[i] / 2))
            return total_moment / total_area
        else:
            return f"ERROR: cannot calculate area because there are {self.widths.shape[0]} width dimensions, and {self.heights.shape[0]} height dimensions."
    
    @property
    def total_height(self):
        return sum(self.heights)

    @property
    def max_width(self):
        return max(self.widths)

    @property
    def gross_inertia(self):
        areas = np.zeros(self.heights.shape[0])
        centers = np.zeros(self.heights.shape[0])
        center_offsets = np.zeros(self.heights.shape[0])
        self_inertias = np.zeros(self.heights.shape[0])
        total_centroid = self.gross_centroid
        running_height = 0
        parallel_axis_terms = np.zeros(self.heights.shape[0])
        for i in range(self.heights.shape[0]):
            self_inertias[i] = self.widths[i] / 12 * self.heights[i] ** 3
            areas[i] = self.widths[i] * self.heights[i]
            centers[i] = self.heights[i] / 2 + running_height
            center_offsets[i] = centers[i] - total_centroid
            running_height += self.heights[i]
            parallel_axis_terms[i] = areas[i] * center_offsets[i] ** 2
        return sum(self_inertias) + sum(parallel_axis_terms)

    @property
    def radius_of_gyration(self):
        return math.sqrt(self.gross_inertia / self.gross_area)
    
    def get_region_count(self, widths_and_heights):
        widths_and_heights = np.array(widths_and_heights)
        return widths_and_heights.shape[0]

    def append_region(self, additional_regions):
        new_cross_section = np.append(self._cross_section, additional_regions)
        new_cross_section.resize(int(new_cross_section.size / 2), 2)
        self._cross_section = new_cross_section
        self._region_count = self.get_region_count(self._cross_section)

    def insert_region(self, index, additional_regions):
        if (index == -1) | (index > self._region_count):
            self._cross_section = np.vstack((self._cross_section, np.array(additional_regions)))
        else:
            new_cross_section = np.insert(self._cross_section, index, additional_regions, axis=0)
            self._cross_section = new_cross_section
        self._region_count = self.get_region_count(self._cross_section)

    def delete_region(self, index):
        if index > self._region_count:
            index = self._region_count
        elif index == 0:
            index = 1
        if self._region_count != 1:
            self._cross_section = np.delete(self._cross_section, index - 1, 0)
            self._region_count = self.get_region_count(self._cross_section)
        
    # Cannot make this a property and make net_distance optional
    def area(self, net_distance = math.inf):
        if self.is_valid:
            net_heights = self.get_net_heights(self.heights, net_distance)
            count = net_heights.shape[0]
            net_area = 0
            if count != 0:
                for i in range(net_heights.shape[0]):
                    net_area += self.widths[i] * net_heights[i]
            return net_area
        else:
            return f"ERROR: cannot calculate area because there are {self.widths.shape[0]} width dimensions, and {self.heights.shape[0]} height dimensions."

    # Cannot make this a property and make net_distance optional
    def centroid(self, net_distance = math.inf):
        if self.is_valid:   
            net_heights = self.get_net_heights(self.heights, net_distance)
            net_area = self.area(net_distance)
            net_moment = 0
            if net_area != 0:
                for i in range(net_heights.shape[0]):
                    net_moment += (self.widths[i] * net_heights[i] * (self.sum_cumulative(net_heights, False)[i] + net_heights[i] / 2))
            return net_moment / net_area
        else:
            return f"ERROR: cannot calculate centroid because there are {self._widths.shape[0]} width dimensions, and {self._heights.shape[0]} height dimensions."

    def get_width_and_height_arrays(self, cross_section):
        cross_section = np.array(cross_section)
        widths = cross_section[:,0]
        heights = cross_section[:,1]
        return widths, heights
        
    def check_if_valid(self, widths, heights):
        return True if widths.shape[0] == heights.shape[0] else False
        
    def get_net_heights(heights, net_distance = math.inf):
        lower_bound_heights = self.sum_cumulative(heights, False)
        net_heights = np.zeros(heights.shape[0])
        for i in range(net_heights.shape[0]):
            net_heights[i] = min(max(0, net_distance - lower_bound_heights[i]), heights[i])
        return net_heights
    
    def sum_cumulative(my_data, include_current: bool = True):
        sum_array = np.zeros(my_data.shape[0])
        if include_current:
            for i in range(sum_array.shape[0]):
                sum_array[i] = np.sum(my_data[0:i + 1])
        else:
            for i in range(1, sum_array.shape[0]):
                sum_array[i] = np.sum(my_data[0:i])     
        return sum_array  
    
    
print(f'{__name__} <version {__version__}> successfully imported')