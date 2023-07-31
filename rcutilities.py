"""Library of general functions for cleaning and preparing data."""
__version__ = "0.0.1"
__author__ = "Ben Fisher"
#__all__ = ['get_list_dimensions', 'is_list', 'to_numpy', 'rescale_value']


import numpy as np

def get_list_dimensions(my_list):
    """Use numpy to get the shape of a python list"""
    if type(my_list) == list:
        temp = np.array(my_list)
        return temp.shape
    else:
        return 'Not a list type'

def is_list(my_list):
    """Check if a variable is a python list type."""
    return True if type(my_list) == list else False

def to_numpy(my_0darray):
    """Creates a numpy array. If 0D array or single value, turns it into a 1D array. Existing numpy arrays are passed through. Python lists are turned into same dimension numpy arrays."""
    if type(my_0darray) == np.ndarray and my_0darray.shape == ():
        return np.array([my_0darray])
    elif type(my_0darray) == int or type(my_0darray) == float:
        return np.array([my_0darray])
    else:
        return np.array(my_0darray)
    
def rescale_value(value, range_min, range_max,
                  new_scale_min = 0.0, 
                  new_scale_max = 1.0, 
                  is_strictly_bounded: bool = True):
    """Rescales a value compared to an original range, and rescales to a new 'scale' range. If strict bounds are imposed, will return the min or max scale limit if the rescaled value is outside the new bounds"""
    rescaled_value = (value - range_min) / (range_max - range_min) * (new_scale_max - new_scale_min) + new_scale_min 
    if is_strictly_bounded == True:
        rescaled_value = min(new_scale_max, max(new_scale_min, rescaled_value))
    return rescaled_value

def get_max_dim(values):
    # TODO: This function needs some work.
    arr = to_numpy(values)
    return len(arr)

######## Below here is a work in progress ########

def region_count(concrete_regions):
    return to_numpy(concrete_regions).shape[0]

def get_region_heights(concrete_regions):
    concrete_regions = to_numpy(concrete_regions)
    count = region_count(concrete_regions)
    heights = np.zeros(count)
    running_total = 0
    for i in range(count):
        running_total += concrete_regions[i][1]
        heights[i] = running_total
    return heights

def get_max_height(concrete_regions):
    return max(get_region_heights(concrete_regions))

def get_lower_bound_region_heights(concrete_regions):
    concrete_regions = to_numpy(concrete_regions)
    count = region_count(concrete_regions)
    heights = np.zeros(count)
    running_total = 0
    for i in range(1, count):
        running_total += concrete_regions[i-1][1]
        heights[i] = running_total
    return heights

def get_region(c, heights):
    region = 0
    for i in range(heights.shape[0]):
        if c <= heights[i]:
            region = i
            break
        elif c >= max(heights):
            region = heights.shape[0] - 1
    return region

def get_distance_into_region(c, concrete_regions):
    region = get_region(c, get_region_heights(concrete_regions))
    lower_bound = get_lower_bound_region_heights(concrete_regions)[region]
    return min(max(0, c - lower_bound), get_max_height(concrete_regions))

def region_areas(concrete_regions, c = math.inf, get_total: bool = False):
    concrete_regions = to_numpy(concrete_regions)
    count = region_count(concrete_regions)
    areas = np.zeros(count)
    for i in range(0, count):
        areas[i] = concrete_regions[i][0] * concrete_regions[i][1]
    return areas if get_total == False else sum(areas)

def get_region_height_bounds(concrete_regions):
    concrete_regions = to_numpy(concrete_regions)
    count = region_count(concrete_regions)
    heights = np.zeros((count, 2))
    running_total = 0
    for i in range(count):
        heights[i,0] = running_total
        running_total += concrete_regions[i][1]
        heights[i,1] = running_total
    return heights

def region_areas_net(concrete_regions, c = math.inf, get_total: bool = False):
    concrete_regions = to_numpy(concrete_regions)
    region_bounds = get_region_height_bounds(concrete_regions)
    count = concrete_regions.shape[0]
    areas = np.zeros(count)
    for i in range(0, count):
        areas[i] = concrete_regions[i][0] * concrete_regions[i][1] * rescale_value(c, region_bounds[i][0], region_bounds[i][1])
    return areas if get_total == False else sum(areas)

def get_region_thicknesses(concrete_regions):
    concrete_regions = to_numpy(concrete_regions)
    return concrete_regions[:, 1]

def get_distance_to_net_centroid(concrete_regions, c = math.inf):
    concrete_regions = to_numpy(concrete_regions)
    thicknesses = get_region_thicknesses(concrete_regions)
    areas = region_areas_net(concrete_regions, c, get_total = False)
    lower_bounds = get_lower_bound_region_heights(concrete_regions)
    moment = 0
    area = 0
    for i in range(thicknesses.shape[0]):
        dist_to_cg = lower_bounds[i] + thicknesses[i] * 0.5
        area += areas[i]
        moment += dist_to_cg * areas[i]
    if area == 0:
        return 0
    else:
        return moment / area


print(f'{__name__} <version {__version__}> successfully imported')