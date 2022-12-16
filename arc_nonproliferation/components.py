from abc import ABC, abstractmethod
import openmc

class Component(ABC):
    """Implement common interface for components"""
    def __init__(self, exclude=True):
        self.exclude = exclude

    def __and__(self, other):
        return openmc.Intersection((self._hull_reg, other))

    def __or__(self, other):
        return openmc.Union((self._hull_reg, other))

    def __invert__(self):
        return ~self._hull_reg

class Blanket(Component):
    def __init__(self, inner_points, outer_points, exclude=True):
        super().__init__()
        self._inner_polygon = openmc.model.Polygon(inner_points,basis='rz')
        self._outer_polygon = openmc.model.Polygon(outer_points,basis='rz')
        self._hull_reg = -self.outer_polygon
        self._body = openmc.Cell(region=(+self._inner_polygon & -self._outer_polygon))
        self._cells = [self._body]
        self._cells.append(openmc.Cell(region=-self._inner_polygon))
        self.exclude = exclude
        
    @property
    def body(self):
        return self._body

    @property
    def cells(self):
        return self._cells

    @property
    def outer_polygon(self):
        return self._outer_polygon

class PFC(Component):
    def __init__(self, thickness, vessel=None, offset=0, lcfs=None, points=None):
        '''
        Defines plasma facing component optionally as:
            1) A first wall of [thickness] offset [offset] from [vessel] inner
               wall. NOTE: Vessel interior cell replaced by PFC interior & gap.
            2) A first wall of [thickness] offset [offset] outside of [lcfs],
               defined as a set of points
            3) A set of [points] defining the poloidal cross-section of the
               component, which will be revolved toroidally. If thickness = 0,
               polygon made from [points] will be filled.
        '''
        
        super.__init__()
        
        # Check inputs
        if vessel is None and points is None and lcfs is None:
            raise ValueError('Neither vessel nor lcfs nor set of points '
                             'specified. Cannot create PFC component.')
        elif vessel is not None:
            if hasattr(vessel, 'inner_polygon'):
                self._outer_polygon = vessel.inner_polygon.offset(-offset)
                self._inner_polygon = self.outer_polygon.offset(-thickness)
                self._area = +self.inner_polygon & -self.outer_polygon
                self._outside = -vessel.inner_polygon & +self.outer_polygon
                self._inside = -self.inner_polygon
                self._cells = [openmc.Cell(region=self.area),
                               openmc.Cell(region=self.inside),
                               openmc.Cell(region=self.outside)]
            else:
                raise TypeError('No inner_polygon attribute defined for '
                                'vessel component passed to PFC component')
        elif lcfs is not None:
            lcfs = openmc.Polygon(lcfs, basis='rz')
            self._inner_polygon = lcfs.offset(offset)
            self._outer_polygon = self.inner_polygon.offset(thickness)
            self._area = +self.inner_polygon & -self.outer_polygon
            self._inside = +lcfs & -self.inner_polygon
            self._outside = None
            self._cells = [openmc.Cell(region=self.area),
                           openmc.Cell(region=self.inside)]
        else:
            base_polygon = openmc.Polygon(points, basis='rz')
            self.area = -base_polygon

            if thickness < 0:
                self._outer_polygon = openmc.Polygon(points, basis='rz').offset(offset)
                self._inner_polygon = self.outer_polygon.offset(thickness)
                self._area = +self.inner_polygon & -self.outer_polygon
                self._inside = -self.inner_polygon
                self._outside = None
                self._cells = [openmc.Cell(region=self.area),
                               openmc.Cell(region=self.inside)]
            elif thickness > 0:
                self._inner_polygon = openmc.Polygon(points, basis='rz').offset(offset)
                self._outer_polygon = self.inner_polygon.offset(thickness)
                self._area = +self.inner_polygon & -self.outer_polygon
                self._inside = -self.inner_polygon
                self._outside = None
                self._cells = [openmc.Cell(region=self.area),
                               openmc.Cell(region=self.inside)]
            else:
                self._outer_polygon = openmc.Polygon(points, basis='rz').offset(offset)
                self._inner_polygon = None
                self._area = -self.outer_polygon
                self._inside = None
                self._outside = None
                self._cells = [openmc.Cell(region=self.area)]
            
        self._cells = [openmc.Cell(region=self.area)]

    @property
    def outer_polygon(self):
        return self._outer_polygon
    
    @property
    def inner_polygon(self):
        return self._inner_polygon
    
    @property
    def area(self):
        return self._area
    
    @property
    def inside(self):
        return self._inside
    
    @property
    def outside(self):
        return self._outside
    
    @property
    def cells(self):
        return self._cells

class Divertor(Component):
    def __init__(self, points, thickness):
        base_polygon = openmc.model.Polygon(points, basis='rz')
        self.inner_polygon = base_polygon
        self.outer_polygon = base_polygon.offset(thickness)
        self._hull_reg = +self.inner_polygon & -self.outer_polygon
        self.inside_reg = -self.inner_polygon
        self.outside_reg = +self.outer_polygon
        self._cells = [openmc.Cell(region=self.inside_reg)]
        self._cells.append(openmc.Cell(region=self._hull_reg))
        self._cells.append(openmc.Cell(region=self.outside_reg))

    @property
    def cells(self):
        return self._cells