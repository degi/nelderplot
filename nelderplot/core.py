"""
Nelder plot designer core module.
The module contains the core classes use for designing the plot.
This module can be used without GUI interfaces.
  
@author: Degi Harja Asmara
"""
from math import log, tan, sqrt, radians, exp, degrees, atan, sin, cos
import simplekml
import utm

class Point:
    """ a class represent the 2D point with x and y coordinates """
    def __init__(self, x = 0.0, y = 0.0):
        self.x = float(x)
        self.y = float(y)
    
    def __add__(self, other):
        if isinstance(other, tuple):
            return Point(self.x + other[0], self.y + other[1])
        if isinstance(other, Point):
            return Point(self.x + other.x, self.y + other.y)
    
    def get(self):
        return self.x, self.y
    
    def __str__(self):
        return str((self.x, self.y)) 

class Spoke:
    """ a class represent the spoke line in Nelder wheel plot """
    epsilon = 0.000001
    
    def __init__(self, nelder_plot, intersect_point: Point = Point(), angle_degree = 0.0, 
                 index = 0, end_point: Point = Point()):
        """ initiate the spoke by their intersected point and angle """
        self.intersect_point = intersect_point
        self.angle_degree = angle_degree
        self.slope = None
        self.intercept = None
        # special case for angle with infinite tangent (parallel  with y axis)
        if not(abs(angle_degree - 90) < self.epsilon or abs(angle_degree - 270) < self.epsilon):
            self.slope = tan(radians(angle_degree))
            self.intercept = intersect_point.y - self.slope * intersect_point.x
        self.end_point = end_point
        self.index = index
        self.nelder_plot = nelder_plot
        self.dataPoints = []
        self.is_border = False
        
    def updateAngle(self, end_point: Point):
        """ update the vector angle relative to x axis defined by and an end point coordinate"""
        self.end_point = end_point
        dy = end_point.y - self.intersect_point.y
        dx = end_point.x - self.intersect_point.x
        a = 0
        if dx < 0: 
            a = 180
        elif dy < 0:
            a = 360
        if dx == 0:
            self.angle_degree = 90 if dy > 0 else 270
        else:    
            self.angle_degree = a + degrees(atan(dy/dx)) 
       
class Wheel:
    """  a class represent the wheel on Nelder plot"""
    epsilon = 0.000001
    
    def __init__(self, nelder_plot, center_point: Point = Point(), radius = 0.0, index = 0):
        """ initiate the wheel by center coordinate and radius """
        self.center_point = center_point
        self.radius = radius
        self.index = index
        self.nelder_plot = nelder_plot
        self.is_border = False
    
    def intersectionCircleLine(self, xc, yc, r, m, i):
        """ 
        calculate the intersection of circle with center coordinates (xc, yc), radius (r)
        and a line with slope (m) and intercept (i) 
        """
        a = m*m + 1
        b = 2 * (m * i - m * yc - xc)
        c = yc*yc - r*r + xc*xc - 2*i*yc + i*i
        sq = b*b - 4*a*c
        if sq < 0: return None
        x1 = (-b + sqrt(sq))/(2*a)
        y1 = m*x1 + i
        x2 = (-b - sqrt(sq))/(2*a)
        y2 = m*x2 + i
        return x1, y1, x2, y2
    
    def intersectionCircleAxis(self, xc, yc, r, x = None, y = None):
        """ 
        calculate the intersection of circle with center coordinates (xc, yc), radius (r)
        and a line parallel with axis line 
        """
        if y is not None:
            a = (r - yc + y)*(r + yc - y)
            if a < 0: return None, None
            if a == 0: return xc, None
            x1 = xc + sqrt(a)
            x2 = xc - sqrt(a)
            return x1, x2
        if x is not None:
            b = (r - xc + x)*(r + xc - x)
            if b < 0: return None, None
            if b == 0: return yc, None
            y1 = yc + sqrt(b)
            y2 = yc - sqrt(b)
            return y1, y2
        return None
    
    def getSpokeIntersection(self, spoke: Spoke):
        """ return the intersection point between spoke line and this wheel """
        if spoke.slope == None:
            if spoke.angle_degree < 180:
                return DataPoint(self.nelder_plot, self, spoke, 
                                 self.center_point.x, self.center_point.y + self.radius)
            else:
                return DataPoint(self.nelder_plot, self, spoke, 
                                 self.center_point.x, self.center_point.y - self.radius)
        if abs(spoke.angle_degree) < self.epsilon:
                return DataPoint(self.nelder_plot, self, spoke, 
                                 self.center_point.x + self.radius, self.center_point.y)
        if abs(spoke.angle_degree - 180) < self.epsilon:
                return DataPoint(self.nelder_plot, self, spoke, 
                                 self.center_point.x - self.radius, self.center_point.y)
        x1, y1, x2, y2 = self.intersectionCircleLine(*self.center_point.get(), self.radius, 
                                                     spoke.slope, spoke.intercept)
        if spoke.angle_degree < 180:
            if y1 >= y2:
                return DataPoint(self.nelder_plot, self, spoke, x1, y1)
            else:
                return DataPoint(self.nelder_plot, self, spoke, x2, y2)  
        else:
            if y1 < y2:
                return DataPoint(self.nelder_plot, self, spoke, x1, y1)
            else:
                return DataPoint(self.nelder_plot, self, spoke, x2, y2)      
    
    def isInsideRange(self, val, maxVal):
        """ check for val validation """
        return not (val == None or val < 0 or val > maxVal)
    
    def getPlotIntersection(self, width, height):
        """ return the intersection of this wheel with the plot rectangle """
        t = ()
        ys = self.intersectionCircleAxis(*self.center_point.get(), self.radius, x = 0)
        if ys[0] is not None and ys[1] is not None:  
            for y in ys:
                if self.isInsideRange(y, height): t += (Point(0, y),)
        ys = self.intersectionCircleAxis(*self.center_point.get(), self.radius, x = width)
        if ys[0] is not None and ys[1] is not None: 
            for y in ys:
                if self.isInsideRange(y, height): t += (Point(width, y),)
        xs = self.intersectionCircleAxis(*self.center_point.get(), self.radius, y = 0)
        if xs[0] is not None and xs[1] is not None: 
            for x in xs:
                if self.isInsideRange(x, width): t += (Point(x, 0),)
        xs = self.intersectionCircleAxis(*self.center_point.get(), self.radius, y = height)
        if xs[0] is not None and xs[1] is not None: 
            for x in xs:
                if self.isInsideRange(x, width): t += (Point(x, height),)
        return t
    
    def getRadius(self):
        """ return the radius property """
        return self.radius
    
    def getSpacing(self):
        """ return the spacing variable """
        return self.nelder_plot.countSpacing(self.getRadius())    
    
    def getDensity(self):
        """ return the density variable """
        return self.nelder_plot.countDensity(self.getRadius())
    
class DataPoint(Point):
    """ an extend of point class with extra information of tree location properties and variables """
    def __init__(self, nelder_plot, wheel: Wheel, spoke: Spoke, x, y, group_id = 0):
        super().__init__(x, y) 
        self.wheel = wheel 
        self.spoke = spoke 
        self.nelder_plot = nelder_plot
        self.group_id = group_id
        self.is_border = False
    
    def getRadius(self):
        """ return the radius variable """
        return self.wheel.getRadius()
    
    def getLabel(self):
        """ return the identification label """
        return "W"+str(self.wheel.index)+"S"+str(self.spoke.index)
    
    def getDescription(self):
        """ return the description for KML pointer information """
        if self.is_border:
            return "Border plant\nSpeciesID = {}\nRadius = {:.5g}".format(
                self.group_id + 1, self.getRadius())
        return "SpeciesID = {}\nSpacing = {:.5g}\nDensity = {:.5g}\nRadius = {:.5g}".format(
            self.group_id + 1, self.getSpacing(), self.getDensity(), self.getRadius()) 
    
    def getSpacing(self):
        """ return the spacing variable """
        return self.wheel.getSpacing()
    
    def getDensity(self):
        """ return the density variable """
        return self.wheel.getDensity()
    
    def __str__(self):
        return str((self.x, self.y)) 

class NelderPlot:
    """ 
    The main class of Nelder plot designer. 
    The class has either standard and practical input parameters. 
    and it will generate Nelder variable including the position coordinates for tree plantation.    
    """
    
    center_config = ("center", "corner", "center_x", "center_y", "custom")
    summary_info = ("Area size (ha)", "Number of plants", "Number of replications")
    param_type_label = ("Practical Parameters", "Basic Parameters", "Center of Rotation", "Number of Species","Summary")
    param_basic_label = ("Plot width (m)", "Plot height (m)", 'Initial radius (m)', 'Rate of change', "Spokes angle (degree)")
    param_practical_label = ("Minimum Density (trees/ha)", "Maximum Density (trees/ha)", "Number of density Ranges", "Rectangularity")

    variable_table_header = ("ID","Radius", "Spacing (m2)", "Density (tress/ha)")
    datapoint_table_header = ("ID", "X", "Y", "Radius", "Spacing (m2)", "Density (tress/ha)", "SpeciesID")
    epsilon_plot = 0.001
    area_unit = 10000
    shapes = ("https://storage.googleapis.com/support-kms-prod/SNP_2752125_en_v0",
              "https://storage.googleapis.com/support-kms-prod/SNP_2752068_en_v0",
              "https://storage.googleapis.com/support-kms-prod/SNP_2752129_en_v0")
    
    def __init__(self, width = 0.0, height = 0.0, init_radius = 0.0, rate_of_change = 0.0, div_angle_degree = 0.0):
        """ initiate the Nelder by standard parameters"""
        self.spokes = []
        self.wheels = {}
        self.datapoints = []
        self.init_radius = 0
        self.number_of_wheels = 0
        self.rate_of_change = 0 
        self.div_angle_degree = 0
        self.non_centrality = 0
        self.rectangularity = 1
        self.center_point =  None
        self.number_of_species = 1
        self.is_alternate_spokes = False
        self.setBasicParameters(width, height, init_radius, rate_of_change, div_angle_degree)
    
    def setDimension(self, width = 0.0, height = 0.0):
        """ define the plot dimension """
        self.width = float(width)
        self.height = float(height)
    
    def setRotationCenter(self, center_point: Point = None, opt = "custom"):
        """ define the center of rotation coordinates directly or using the pre-defined options (opt) """
        if opt == self.center_config[4] and not center_point == None:
            self.center_point = center_point
        elif opt == self.center_config[1]:
            self.center_point = Point(0, 0)
        elif opt == self.center_config[2]:
            self.center_point = Point(self.width/2, 0)
        elif opt == self.center_config[3]:
            self.center_point = Point(0, self.height/2)
        else:
            self.center_point = Point(self.width/2, self.height/2)    
            
    def setBasicParameters(self, width = 0.0, height = 0.0, init_radius = 0.0, rate_of_change = 0.0, 
                           div_angle_degree = 0.0, opt_center = center_config[0], center_point: Point = None):
        """ set the Nelder by set of standard parameters """
        self.setDimension(width, height)
        self.setRotationCenter(center_point, opt_center)
        self.init_radius = float(init_radius)
        self.rate_of_change = float(rate_of_change)
        self.div_angle_degree = float(div_angle_degree)
        self.updateNonCentrality(self.rate_of_change)
        if self.init_radius > 0 and self.rate_of_change > 0:
            rMax = max(*self.center_point.get(), 
               abs(self.width - self.center_point.x), abs(self.height - self.center_point.y))
            self.number_of_wheels = int(log(rMax/self.init_radius)/log(self.rate_of_change))
            self.rectangularity = radians(self.div_angle_degree)/ (sqrt(self.rate_of_change) - 1/sqrt(self.rate_of_change))
        self.updatePlotProperties()
    
    def setPracticalParameters(self, min_density, max_density, n_range, rectangularity = 1, 
                               opt_center = center_config[0], center_point: Point = None):
        """ set the Nelder by set of practical parameters """
        if not(min_density > 0 and max_density > 0 and n_range > 0): return 
        self.number_of_wheels = int(n_range) + 1
        if n_range == 1: return
        a = exp((log(max_density) - log(min_density))/(2 * n_range - 2))
        self.rate_of_change = a
        self.updateNonCentrality(a)
        if rectangularity == 0: rectangularity = 1
        self.rectangularity = rectangularity
        t = rectangularity * (sqrt(a) - 1/sqrt(a))
        self.div_angle_degree = degrees(t)
        self.init_radius = sqrt(20000/(max_density * t * (a**3 - a)))
        r = self.countRadius(n_range + 1)
 
        if opt_center == self.center_config[0]:
            self.setDimension(r*2, r*2)   
        elif opt_center == self.center_config[1]: 
            self.setDimension(r, r)
        elif opt_center == self.center_config[2]: 
            self.setDimension(r*2, r)
        elif opt_center == self.center_config[3]: 
            self.setDimension(r, r*2)
        elif opt_center == self.center_config[4] and not center_point == None:
            self.setDimension(center_point.x + r, center_point.y + r)  
        self.setRotationCenter(center_point, opt_center)
        
        self.updatePlotProperties()
    
    def getBasicParameters(self):
        """ return the standard parameters """
        return self.width, self.height, self.init_radius, self.rate_of_change, self.div_angle_degree
    
    def getPracticalParameters(self):
        """ return the practical parameters """
        return self.getMinimumDensity(), self.getMaximumDensity(), self.number_of_wheels-1, self.rectangularity
        
    def updateNonCentrality(self, a):
        """ update the non centrality parameter """
        if a > 0 and a is not 1:
            self.non_centrality = 100*(2 - (a + 1/a))/(2*(a - 1/a))

    def isValid(self):
        """ check the parameter validity """
        if self.init_radius == 0 or self.rate_of_change == 0: return False
        return True
    
    def clearData(self):
        """ reset the plot design and data """
        self.spokes.clear()
        self.wheels.clear()
        self.datapoints.clear()
        
    def updatePlotProperties(self):
        """ 
        compute the Nelder plot properties and variables.
        this function is called automatically after setting the parameters.
        """
        self.clearData()
        if not self.isValid(): return
        self.number_of_species = max(min(self.number_of_species, 3), 1)
        # create the wheels
        for i in range(self.number_of_wheels + 1):
            r = self.countRadius(i)  
            self.wheels[i] = Wheel(self, self.center_point, r, i)
            self.wheels[i].group_id = i % self.number_of_species 
        if self.div_angle_degree == 0: return
        # flag the innermost and outermost wheels as border sampling
        self.wheels[0].is_border = True
        self.wheels[self.number_of_wheels].is_border = True
        # check for wheel intersection with the plot
        pList = None
        if self.isInside(self.center_point):
            pList = self.wheels[self.number_of_wheels].getPlotIntersection(self.width, self.height)
        else:
            pList = self.wheels[0].getPlotIntersection(self.width, self.height)
        spokeList = ()
        for p in pList:
            spoke = Spoke(self, self.center_point)
            spoke.updateAngle(p)
            spokeList += (spoke,)
        # define the start angle based on wheel intersection with the plot
        startAngles = ()
        for sp in spokeList:
            # check the validity of starting angle by creating the spoke with addition of division angle
            # if the spoke inside the plot, then it will be a valid starting angle
            spoke = Spoke(self, self.center_point, (sp.angle_degree + self.div_angle_degree) % 360, i)
            p0 = self.wheels[0].getSpokeIntersection(spoke)
            pn = self.wheels[self.number_of_wheels].getSpokeIntersection(spoke)
            if self.isInside(pn) and self.isInside(p0):
                startAngles += (sp.angle_degree,)
        # create the spokes
        nAngle, modAngle = divmod(360, self.div_angle_degree)
        isBorderSpokes = False if modAngle == 0 and len(startAngles) == 0 else True
        nAngle = int(nAngle)
        startAngle = 0
        iStart = 0
        if len(startAngles) > 0: startAngle = startAngles[iStart]
        a = 0
        for i in range(nAngle + 1):
            angleDeg = (startAngle + self.div_angle_degree * a) % 360
            spoke = Spoke(self, self.center_point, angleDeg, i)
            spoke.is_border = isBorderSpokes and (a == 0 or i == nAngle) 
            # get the intersection of spoke line and wheel
            p0 = self.wheels[0].getSpokeIntersection(spoke)
            pn = self.wheels[self.number_of_wheels].getSpokeIntersection(spoke)
            a += 1
            if self.isInside(p0) and self.isInside(pn):
                # if the intersection point is inside the plot, put it on the containers
                spoke.end_point = pn
                spoke.group_id = (i % self.number_of_species)
                self.spokes.append(spoke)
                p0.is_border = pn.is_border = True 
                p0.group_id = self.countGroupID(self.wheels[0], spoke)
                pn.group_id = self.countGroupID(self.wheels[self.number_of_wheels], spoke)
                self.addPoint(pn)
                self.addPoint(p0)
                # calculate for spoke intersection with the other wheels
                for j in range(1, self.number_of_wheels):
                    pj = self.wheels[j].getSpokeIntersection(spoke)
                    pj.is_border = spoke.is_border
                    pj.group_id = self.countGroupID(self.wheels[j], spoke)
                    self.addPoint(pj)
            else:
                # if the recent intersection is outside of the plot and it is the later checking
                # then previous spoke should be a border spoke  
                if len(self.spokes) > 1:
                    self.spokes[len(self.spokes) - 1].is_border = True
                    for d in self.spokes[len(self.spokes) - 1].dataPoints: d.is_border = True
                # reset the start angle if it has more then one start angle list
                if len(startAngles) > 1 and iStart < 1: 
                    iStart += 1
                    startAngle = startAngles[iStart]
                    i -= 1
                    a = 0
                else:
                    break
    
    def countGroupID(self, w: Wheel, s: Spoke):
        """ compute group id or species id for alternate spokes and/or wheels """
        if self.is_alternate_spokes: return s.group_id
        return (w.group_id + s.group_id) % self.number_of_species
            
    def countRadius(self, index):
        """ count the radius by the wheel index """
        if index < 0: return 0
        return self.init_radius * self.rate_of_change ** index
   
    def addPoint(self, point:DataPoint):
        """ add the data point of tree position to the container """
        self.datapoints.append(point)
        point.spoke.dataPoints.append(point)
    
    def countSpacing(self, r):
        """ count the spacing by the defined radius """
        if self.rate_of_change == 0: return 0
        return r**2 * radians(self.div_angle_degree) * (self.rate_of_change - 1/self.rate_of_change)/2
                                
    def isInside(self, point: Point):
        """ check the position whether it is inside the plot"""
        if (point.x >= -self.epsilon_plot and point.y >= -self.epsilon_plot 
            and point.x <= self.width + self.epsilon_plot and point.y <= self.height + self.epsilon_plot): 
            return True
        return False
    
    def getMinimumSpacing(self):
        """ return minimum spacing """
        return self.countSpacing(self.countRadius(1))
    
    def getMaximumSpacing(self):
        """ return maximum spacing """
        return self.countSpacing(self.countRadius(self.number_of_wheels-1))
    
    def getMaximumDensity(self):
        """ return maximum plot density """
        s = self.getMinimumSpacing() 
        if s <= 0: return 0
        return self.area_unit/s

    def getMinimumDensity(self):
        """ return minimum plot density """
        s = self.getMaximumSpacing()
        if s <= 0: return 0
        return self.area_unit/s
    
    def countDensity(self, r):
        """ calculate the plot density """
        s = self.countSpacing(r)
        if s <= 0: return 0
        return self.area_unit/s
    
    def getReplications(self):
        """ calculate the experiment unit replication  """
        n = 0
        for sp in self.spokes:
            if not sp.is_border: n += 1
        return n // self.number_of_species
    
    def getSummaryInfo(self):
        """ return the summary of defined plot design """
        summary = {}
        summary[self.summary_info[0]] = (self.width * self.height)/self.area_unit 
        summary[self.summary_info[1]] = len(self.datapoints)
        summary[self.summary_info[2]] = self.getReplications()
        return summary
    
    def getDescription(self):
        """ return the description for KML plot data file """
        desc = ""
        desc += "# {} \n".format(self.param_type_label[0])
        p_pract = self.getPracticalParameters()
        for t, v in zip(self.param_practical_label, p_pract):
            desc += "{} = {:.5g}\n".format(t, v)
        desc += "\n# {} \n".format(self.param_type_label[1])
        p_basic = self.getBasicParameters()
        for t, v in zip(self.param_basic_label, p_basic):
            desc += "{} = {:.5g}\n".format(t, v) 
        desc += "\n# {} \n".format(self.param_type_label[2])
        desc += "({:.5g}, {:.5g})\n".format(self.center_point.x, self.center_point.y)
        desc += "\n# {} \n".format(self.param_type_label[4])
        desc += "{} = {:.5g}\n".format(self.param_type_label[3], self.number_of_species)
        info = self.getSummaryInfo()
        for t in self.summary_info:
            desc += "{} = {:.5g}\n".format(t, info[t])
        desc += "\n# Plot variables\n"
        for h in self.variable_table_header: desc += "{}; ".format(h)
        desc += "\n"
        vt = self.getVariableTable()
        for row in vt: 
            for v in row: 
                if isinstance(v, str):
                    desc += "{}; ".format(v)
                else:
                    desc += "{:.5g}; ".format(v)
            desc += "\n"
        return desc
    
    def getVariableTable(self):
        """ 
        return the data table of plot variables 
        header = "ID","Radius", "Spacing (m2)", "Density (tress/ha)"
        """
        s = []
        f = '{:.5g}'
        for i, w in self.wheels.items():
            sp = "" if w.is_border else f.format(w.getSpacing())
            de = "" if w.is_border else f.format(w.getDensity())
            s.append(("W" + str(i), f.format(w.getRadius()), sp, de))
        return s
            
    def getDatapointsTable(self):
        """ 
        return the data table of position point     
        header = "ID", "X", "Y", "Radius", "Spacing (m2)", "Density (tress/ha)", "SpeciesID"
        """
        s = []
        f = '{:.5g}'
        for p in self.datapoints:
            s.append((p.getLabel(), f.format(p.x), f.format(p.y), f.format(p.getRadius()), 
                      "" if p.is_border else f.format(p.getSpacing()), 
                      "" if p.is_border else f.format(p.getDensity()), p.group_id+1))
        return s
    
    def getPlotPolygon(self):
        """ return the polygon coordinates """
        return [(0,0), (self.width, 0), (self.width, self.height), (0, self.height)]
    
    
    def rotatePoint(self, x, y, rotationDegree):
        """ transforms the coordinates by rotation from origin """
        sinR = sin(radians(rotationDegree))
        cosR = cos(radians(rotationDegree))
        xr = x * cosR + y * sinR
        yr = y * cosR - x * sinR  
        return xr, yr

    def getKML(self, latitude, longitude, rotation_degree):
        """ 
        return the KML format text with plot origin coordinates on (latitude, longitude)
        and rotate by defined degrees
        """
        plot_utm = utm.from_latlon(latitude, longitude)
        kml = simplekml.Kml()
        pol = kml.newpolygon(name='Nelder Plot')
        ps = self.getPlotPolygon()
        b = []
        for p in ps:
            pr = p
            if not rotation_degree == 0: pr = self.rotatePoint(*pr, rotation_degree)
            latlon = utm.to_latlon(plot_utm[0] + pr[0], plot_utm[1] + pr[1], plot_utm[2], plot_utm[3])
            b.append((latlon[1], latlon[0]))
        pol.outerboundaryis = b
        pol.style.polystyle.color = simplekml.Color.changealphaint(10, simplekml.Color.purple)
        pol.style.linestyle.width = 3
        pol.description = self.getDescription()
        for p in self.datapoints:
            pr = (p.x, p.y)
            if not rotation_degree == 0: pr = self.rotatePoint(*pr, rotation_degree)
            latlon = utm.to_latlon(plot_utm[0] + pr[0], plot_utm[1] + pr[1], plot_utm[2], plot_utm[3])
            pkml = kml.newpoint(name=p.getLabel(),  coords=[(latlon[1], latlon[0])])
            pkml.description = p.getDescription()
            if p.is_border:
                pkml.style.iconstyle.icon.href = self.shapes[p.group_id] 
            else:
                pkml.style.iconstyle.icon.href = self.shapes[p.group_id]
            pkml.style.iconstyle.scale = 0.3
        return kml.kml()
