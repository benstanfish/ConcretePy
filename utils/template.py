import math
from typing import (
    List,
    Dict,
    Tuple,
    Literal,
    overload
)

class Point():

    def __init__(self, x: float, y: float) -> None:
        self.x = x
        self.y = y
        
    @property
    def coords(self) -> Tuple[float, float]:
        return (self.x, self.y)


class BoundingBox():

    def __init__(self, *args, **kwargs):
        self.points = []
        if len(args) > 0:
            if len(args) == 4:
                self.points.clear
                self.points.append(Point(args[0], args[1]))
                self.points.append(Point(args[2], args[1]))
                self.points.append(Point(args[2], args[3]))
                self.points.append(Point(args[0], args[3]))
            elif len(args) == 2 and isinstance(args[0], Point):
                self.points.clear
                self.points.append(args[0])
                self.points.append(Point(args[1].x, args[0].y))
                self.points.append(args[1])
                self.points.append(Point(args[0].x, args[1].y))
        if len(kwargs) > 0:
            self.points.clear
            self.points.append(Point(kwargs['left'], kwargs['top']))
            self.points.append(Point(kwargs['right'], kwargs['top']))
            self.points.append(Point(kwargs['right'], kwargs['bottom']))
            self.points.append(Point(kwargs['left'],kwargs['bottom']))

    @property
    def coords(self) -> List[Tuple[float, float]]:
        return [(float(point.x), float(point.y)) for point in self.points]
    
    @property
    def centroid(self) -> Point:
        xs = [point.x for point in self.points]
        ys = [point.y for point in self.points]
        return Point(sum(xs)/len(xs), sum(ys)/len(ys))

    @property
    def area(self) -> float:
        """Area calculated using the Trapezoid variant of the Shoelace Formula."""
        xs = [point.x for point in self.points]
        ys = [point.y for point in self.points]
        temp = 0
        for i in range(len(xs)):
            temp += (ys[i-1] + ys[i]) * (xs[i-1] - xs[i])
        temp *= 0.5
        return temp

    @property
    def left(self) -> float:
        return min([point.x for point in self.points])

    @property
    def right(self) -> float:
        return max([point.x for point in self.points])

    @property
    def top(self) -> float:
        return max([point.y for point in self.points])

    @property
    def bottom(self) -> float:
        return min([point.y for point in self.points])

    def translate(self, **kwargs) -> None:
        if 'left' in kwargs.keys():
            for point in self.points:
                point.x -= float(kwargs['left'])
        if 'right' in kwargs.keys():
            for point in self.points:
                point.x += float(kwargs['right'])
        if 'up' in kwargs.keys():
            for point in self.points:
                point.y -= float(kwargs['up'])
        if 'down' in kwargs.keys():
            for point in self.points:
                point.y += float(kwargs['down'])
        if 'horizontal' in kwargs.keys():
            for point in self.points:
                point.x += float(kwargs['horizontal'])          
        if 'vertical' in kwargs.keys():
            for point in self.points:
                point.y -= float(kwargs['vertical'])     

    def rotate(self, **kwargs) -> None:
        if list(kwargs.keys())[0] == 'degrees':
            radians = float(math.radians(kwargs['degrees']))
        else:
            radians = float(kwargs['radians'])
        temp = [self._point_rotation(point, radians) for point in self.points]
        self.points.clear
        self.points = temp

    def scale(self, **kwargs) -> None:
        scale = 1
        scale_x = scale
        scale_y = scale
        if 'scale' in kwargs.keys():
            scale = kwargs['scale']
            scale_x = scale
            scale_y = scale
        if 'scale_x' in kwargs.keys():
            scale_x=kwargs['scale_x']
        if 'scale_y' in kwargs.keys():
            scale_y=kwargs['scale_y']
        temp = [self._point_scale(point, scale_x, scale_y) for point in self.points] 
        self.points.clear
        self.points = temp

    def _point_rotation(self, point: Point, radians: float) -> Point:
        """Uniform rotation transformation about the bounding boxes centroid."""
        a, b = self.centroid.coords
        x, y = point.coords
        return Point((x - a) * math.cos(radians) - (y - b) * math.sin(radians) + a, 
                     (x - a) * math.sin(radians) + (y - b) * math.cos(radians) + b)

    def _point_scale(self, point: Point, scale_x: float, scale_y: float) -> Point:
        """Uniform scale transformation about the bounding boxes centroid."""
        a, b = self.centroid.coords
        x, y = point.coords
        return Point((x - a) * scale_x + a, (y - b) * scale_y + b)




class BoundingPolygon():

    def __init__(self, *args):
        self.points = []
        if isinstance(args[0], Tuple):
            self.points.clear
            for coords in args:
                self.points.append(Point(coords[0], coords[1]))
        if isinstance(args[0], Point):
            self.points.clear
            for point in args:
                self.points.append(point)

    @property
    def coords(self) -> List[Tuple[float, float]]:
        return [(float(point.x), float(point.y)) for point in self.points]
    
    @property
    def mean_center(self) -> Point:
        """Returns the point of the mean center - not the the centroid."""
        xs = [point.x for point in self.points]
        ys = [point.y for point in self.points]
        return Point(sum(xs) / len(xs), sum(ys) / len(ys))

    @property
    def centroid(self) -> Point:
        """Returns the centroid or center of mass, assuming uniform density."""
        a = self.area
        xs = [point.x for point in self.points]
        ys = [point.y for point in self.points]
        x_bar = 0
        y_bar = 0
        for i in range(len(xs)):
            x_bar += (xs[i-1] + xs[i]) * (xs[i - 1] * ys[i] - xs[i] * ys[i - 1])
            y_bar += (ys[i-1] + ys[i]) * (xs[i - 1] * ys[i] - xs[i] * ys[i - 1])
        return Point(x_bar / (6 * a), y_bar / (6 * a))

    @property
    def mbr_center(self) -> Point:
        """Center point of the minimum bounding rectangle."""
        return Point((self.right - self.left) / 2, (self.top - self.bottom) / 2)

    @property
    def area(self) -> float:
        """Area calculated using the Trapezoid variant of the Shoelace Formula."""
        xs = [point.x for point in self.points]
        ys = [point.y for point in self.points]
        temp = 0
        for i in range(len(xs)):
            temp += (ys[i-1] + ys[i]) * (xs[i-1] - xs[i])
        temp *= 0.5
        return abs(temp)

    @property
    def left(self) -> float:
        return min([point.x for point in self.points])

    @property
    def right(self) -> float:
        return max([point.x for point in self.points])

    @property
    def top(self) -> float:
        return max([point.y for point in self.points])

    @property
    def bottom(self) -> float:
        return min([point.y for point in self.points])

    def translate(self, **kwargs) -> None:
        if 'left' in kwargs.keys():
            for point in self.points:
                point.x -= float(kwargs['left'])
        if 'right' in kwargs.keys():
            for point in self.points:
                point.x += float(kwargs['right'])
        if 'up' in kwargs.keys():
            for point in self.points:
                point.y -= float(kwargs['up'])
        if 'down' in kwargs.keys():
            for point in self.points:
                point.y += float(kwargs['down'])
        if 'horizontal' in kwargs.keys():
            for point in self.points:
                point.x += float(kwargs['horizontal'])          
        if 'vertical' in kwargs.keys():
            for point in self.points:
                point.y -= float(kwargs['vertical'])     

    def rotate(self, **kwargs) -> None:
        if list(kwargs.keys())[0] == 'degrees':
            radians = float(math.radians(kwargs['degrees']))
        else:
            radians = float(kwargs['radians'])
        temp = [self._point_rotation(point, radians) for point in self.points]
        self.points.clear
        self.points = temp

    def scale(self, **kwargs) -> None:
        scale = 1
        scale_x = scale
        scale_y = scale
        if 'scale' in kwargs.keys():
            scale = kwargs['scale']
            scale_x = scale
            scale_y = scale
        if 'scale_x' in kwargs.keys():
            scale_x=kwargs['scale_x']
        if 'scale_y' in kwargs.keys():
            scale_y=kwargs['scale_y']
        temp = [self._point_scale(point, scale_x, scale_y) for point in self.points] 
        self.points.clear
        self.points = temp

    def _point_rotation(self, point: Point, radians: float) -> Point:
        """Uniform rotation transformation about the bounding boxes centroid."""
        a, b = self.centroid.coords
        x, y = point.coords
        return Point((x - a) * math.cos(radians) - (y - b) * math.sin(radians) + a, 
                     (x - a) * math.sin(radians) + (y - b) * math.cos(radians) + b)

    def _point_scale(self, point: Point, scale_x: float, scale_y: float) -> Point:
        """Uniform scale transformation about the bounding boxes centroid."""
        a, b = self.centroid.coords
        x, y = point.coords
        return Point((x - a) * scale_x + a, (y - b) * scale_y + b)
    


# bp = BoundingPolygon(Point(3, 0),
#                      Point(14, 3),
#                      Point(10, 4),
#                      Point(13, 11),
#                      Point(4, 13),
#                      Point(0, 8))



bp = BoundingPolygon(Point(0, 0),
                     Point(0, 100),
                     Point(100, 100),
                     Point(100, 0))

print(bp.area)
print(bp.mean_center.coords)
print(bp.centroid.coords)
print(bp.mbr_center.coords)