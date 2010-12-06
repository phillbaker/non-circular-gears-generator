#!/usr/bin/python2.2

from slfmaker import SLFMaker, SLFMakerSphereTeeth
from math import sin,cos,pi,sqrt

class oval (SLFMaker):
    "Makes elliptical gears."
    def __init__(self, teeth, ts, a, e, nodes, iradius, depth, tolerance):
        # various oval parameters
        self.a = a;
        self.nodes = nodes
        
        # make sure eccentricity value is valid
        
        self.c = e * a
        self.e = e
        assert self.e <= 1.0 and self.e >= 0.0
        self.p = a * (1.0 - self.e*self.e)
        self.b = a*a - self.c*self.c

        self.centerOffset = 0.0

        self.iradius = iradius # inner radius, for the hole
        print "Elliptical Gear"
        print "a =", self.a, "c =", self.c, "e =", self.e
        print "teeth =", teeth, "thickness =", depth, "inner radius =", iradius

        # oval values are needed by perimeter, which is needed by init
        SLFMaker.__init__(self, teeth, ts, depth, tolerance)

    def width(self):
        return 2.0*(self.a+self.adendumd+self.dedendumd)

    def outerradius(self,theta):
         return self.p / (1.0 - self.e * cos(self.nodes*theta))

    def innerradius(self, theta):
        return self.iradius

    def dx(self, t):
        return (-self.e*self.nodes*self.p*cos(t)*sin(self.nodes*t))/(1.0-self.e*cos(self.nodes*t))**2 - \
               (self.p*sin(t))/(1.0 - self.e*cos(self.nodes*t))
    
    def dy(self, t):
        return self.p*cos(t)/(1.0-self.e*cos(self.nodes*t)) - \
               self.e*self.nodes*self.p*sin(t)*sin(self.nodes*t) / (1.0-self.e*cos(self.nodes*t))**2

    def dx2(self, t):
        return (2.0*self.e*self.e*self.nodes*self.nodes*self.p*cos(t)*sin(self.nodes*t)**2) / \
               (1.0-self.e*cos(self.nodes*t))**3 + \
               (2.0*self.e*self.nodes*self.p*sin(t)*sin(self.nodes*t)) / \
               (1.0 - self.e*cos(self.nodes*t))**2 - \
               (self.p*cos(t))/(1.0-self.e*cos(self.nodes*t)) - \
               (self.e*self.nodes*self.nodes*self.p*cos(t)*cos(self.nodes*t)) / \
               (1.0-self.e*cos(self.nodes*t))**2

    def dy2(self, t):
        return (2.0*self.e*self.e*self.nodes*self.nodes*self.p*sin(t)*sin(self.nodes*t)**2) / \
               (1.0-self.e*cos(self.nodes*t))**3 - \
               (2.0*self.e*self.nodes*self.p*cos(t)*sin(self.nodes*t)) / \
               (1.0 - self.e*cos(self.nodes*t))**2 - \
               (self.p*sin(t))/(1.0-self.e*cos(self.nodes*t)) - \
               (self.e*self.nodes*self.nodes*self.p*sin(t)*cos(self.nodes*t)) / \
               (1.0-self.e*cos(self.nodes*t))**2        

    def perimeter(self):
        "Computes the perimeter of an oval using Ramanujan's approximation."
        h = ((self.a-self.b)/(self.a+self.b))**2.0
        return pi * (self.a+self.b) * (1.0+(3.0 * h)/(10.0+sqrt(4.0-3.0*h)))

e = oval(30, 10, 1.0, 0.15, 2, 0.0625, 0.25, 0.0001)

e.write('ncgear3.slf')
#for x in range(1,9):
#    print 2.0*pi/x, e.radiusOfCurvature(2.0*pi/x)
