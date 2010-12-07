#!/usr/bin/python

from slfmaker import SLFMaker, SLFMakerSphereTeeth
from math import sin,cos,pi,sqrt

class circle (SLFMaker):
    def __init__(self, teeth, ts, oradius, iradius, depth, tolerance):
        self.centerOffset = 0
        self.iradius = iradius # inner radius, for the hole
        self.oradius = oradius
        print "Circle"
        print "radius =", oradius
        print "teeth =", teeth, "thickness =", depth, "inner radius =", iradius

        # ellipse values are needed by perimeter, which is needed by init
        SLFMaker.__init__(self, teeth, ts, depth, tolerance)

    def outerradius(self,theta):
        # x = a * cos (theta)
        # y = a * sin (theta)
        return self.oradius

    def innerradius(self, theta):
        return self.iradius

    def dx(self, t):
        return -self.oradius * sin(t)

    def dy(self, t):
        return self.oradius * cos(t)

    def dx2(self, t):
        return -self.oradius * cos(t)

    def dy2(self, t):
        return -self.oradius * sin(t)
    
    def perimeter(self):
        "Computes the perimeter of an ellipse using Ramanujan's approximation."
        return self.oradius*pi
        
    def width(self):
        return self.outerradius(0)*2

e = circle(20, 10, 1.0, 0.125, 0.2, 0.001)

#import profile
#profile.run("e.write('ncgear2.slf')")
e.write('cgear4.dxf')
#for x in range(1,9):
#    print 2.0*pi/x, e.radiusOfCurvature(2.0*pi/x)

