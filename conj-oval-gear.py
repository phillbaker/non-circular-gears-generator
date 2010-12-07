#!/usr/bin/python

from conj_slfmaker import ConjugateSLFMaker, ConjugateSLFMakerSphereTeeth
from math import sin,cos,pi,sqrt

class oval (ConjugateSLFMaker):
    "Makes oval gears."
    def __init__(self, teeth, ts, a, e, nodes, period, holedistance, iradius, depth, tolerance):
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
        print "Oval Gear Pair"
        print "a =", self.a, "c =", self.c, "e =", self.e
        print "teeth =", teeth, "thickness =", depth, "inner radius =", iradius

        # oval values are needed by perimeter, which is needed by init
        ConjugateSLFMaker.__init__(self, teeth, ts, period, holedistance, depth, tolerance)

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

    def drdt(self, t):
        return -self.e*self.nodes*self.p*sin(self.nodes*t)/ \
               (1.0 - self.e*cos(self.nodes*t))**2

    def perimeter(self):
        "Computes the perimeter of an oval using Ramanujan's approximation."
        h = ((self.a-self.b)/(self.a+self.b))**2.0
        return pi * (self.a+self.b) * (1.0+(3.0 * h)/(10.0+sqrt(4.0-3.0*h)))

e = oval(40, 10, 1.0, 0.15, 2, 2, 3.0, 0.0625, 0.25, 0.001)
#e = oval(20, 10, 1.0, 0.65, 1, 1, 2.0, 0.0625, 0.25, 0.001)
e.write('ncgear1.dxf')
#for x in range(1,9):
#    print 2.0*pi/x, e.radiusOfCurvature(2.0*pi/x)
