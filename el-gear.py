#!/usr/bin/python

from slfmaker import SLFMaker
from offsetpair_slfmaker import OffsetPairSLFMaker
from math import sin,cos,pi,sqrt

def makeEllipse(teeth, ts, a, e, iradius, depth, tolerance):
    if (teeth % 2 == 0): # even
        cl = OffsetPairSLFMaker
    else:
        cl = SLFMaker

    class ellipse (cl):
        "Makes elliptical gears."
        def __init__(self, teeth, ts, a, e, iradius, depth, tolerance):
            # various ellipse parameters
            self.a = a;
        
            # make sure eccentricity value is valid
        
            self.c = e * a
            self.e = e
            assert self.e <= 1.0 and self.e >= 0.0
            self.p = a * (1.0 - self.e*self.e)
            self.b = a*a - self.c*self.c

            self.iradius = iradius # inner radius, for the hole

            # required variables for centering of the viewing frustum

            self.centerOffset = self.c
        
            print "Elliptical Gear"
            print "a =", self.a, "c =", self.c, "e =", self.e
            print "teeth =", teeth, "thickness =", depth, "inner radius =", iradius

            # ellipse values are needed by perimeter, which is needed by init
            cl.__init__(self, teeth, ts, depth, tolerance)

        def width(self):
            return (self.a + self.adendumd)*2.0

        def outerradius(self,theta):
            return self.p / (1.0 + self.e * cos(theta))

        def innerradius(self, theta):
            return self.iradius

        def dx(self, t):
            return self.e*self.p*cos(t)*sin(t)/(self.e*cos(t)+1.0)**2 - \
                   self.p*sin(t)/(self.e*cos(t)+1.0)

        def dy(self, t):
            return self.e*self.p*sin(t)*sin(t)/(self.e*cos(t)+1.0)**2 + \
                   self.p*cos(t)/(self.e*cos(t)+1.0)
        
        def dx2(self, t):
            return -2.0*self.e*self.p*sin(t)*sin(t) / (self.e*cos(t)+1.0)**2 + \
                   2.0*self.e*self.e*self.p*cos(t)*sin(t)*sin(t)/(self.e*cos(t)+1.0)**3 - \
                   self.p*cos(t)/(self.e*cos(t)+1.0) + \
                   self.e*self.p*cos(t)*cos(t) / (self.e*cos(t)+1.0)**2
        
        def dy2(self, t):
            return 2.0 * self.e*self.e*self.p*sin(t)**3/(self.e*cos(t)+1.0)**3 - \
                   self.p*sin(t) / (self.e*cos(t)+1.0) + \
                   3.0*self.e*self.p*cos(t)*sin(t) / (self.e * cos(t)+1.0)**2
    
        def perimeter(self):
            "Computes the perimeter of an ellipse using Ramanujan's approximation."
            h = ((self.a-self.b)/(self.a+self.b))**2.0
            return pi * (self.a+self.b) * (1.0+(3.0 * h)/(10.0+sqrt(4.0-3.0*h)))

    e = ellipse(teeth, ts, a, e, iradius, depth, tolerance)
    if (teeth % 2 == 0):
        e.secondOffset = a
    return e

e = makeEllipse(20, 10, 1.0, 0.65, 0.0625, 0.25, 0.0001)
e.write('ncgear1.slf')
