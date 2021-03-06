#!/usr/bin/python

from math import sin, cos, pi, sqrt, atan

def involute(a, t):
    """Returns the 2D point [x,y] for the involute of a circle
    with radius a, with the parameter t."""
    return [a * (sin(t) - t*cos(t)), a * (cos(t) + t*sin(t))-a]
#    return [a * (cos(t) + t*sin(t)), a * (sin(t) - t*cos(t))]

def minIndex(l):
    m = l[0]
    mi = 0
    for n in range(1,len(l)):
        if l[n] < m:
            m = l[n]
            mi = n
    return mi

def findClosest(l, v):
    d = map(lambda x: abs(x[0]-v), l)
    return l[minIndex(d)][1]

def findClosestDown(l, v):
    best = filter(lambda x: v - x[0] >= 0.0, l)[0]
    for x in range(len(l)):
        i = v - l[x][0]
        if (i >= 0.0 and i <= best):
            bestIndex = x
            best = i
    return l[bestIndex][1]

def normalize(x, y):
    if (x==0.0 and y==0.0): return (0.0,0.0)
    return (x/sqrt(x*x+y*y),y/sqrt(x*x+y*y))

def toCartesian(r,t):
    return (r*cos(t), r*sin(t))

class SLFMaker:
    def __init__(self, teethCount, ts=5, depth=0.1, tolerance = 0.001):
        """self.perimeter() needs to be defined.
        """
        self.depth = depth
        self.teethLoc = []
        self.teethCount = teethCount
        self.cpitch = self.perimeter() / teethCount
        self.toothSlices = ts
        print "est. circular pitch:", self.cpitch, "est. perimeter: ", self.perimeter()

        self.tolerance = tolerance
        
    def write(self, filename):
        if not self.teethLoc:
            self.calcPoints()
        f = open(filename, 'w')
        self.preamble(f)
        self.amble(f)
        self.postamble(f)        
        f.close()

    def amble(self, f):
        """Takes the pointlist and produces the necessary dxf code for
        the points and the line they're in"""
        import string

        n = 0
        inc = 0.01

        print "teeth completed: ",
        
        teeth_ends = []
        inner_pts = []
        #do this for each tooth:
        for d in self.teethLoc:
            rc = abs(d['rc'])
            inverseInvoluteTableY = []
            inverseInvoluteTableX = []
            t = 0.0
            while t < pi / 2.0:
                i = involute(rc, t)
                inverseInvoluteTableX.append((i[0],t))
                inverseInvoluteTableY.append((i[1],t))
                t += inc

            # t parameter at pitch circle
            td = findClosest(inverseInvoluteTableY,self.dedendumd)

            # t parameter at top of tooth
            tt = findClosest(inverseInvoluteTableY,self.dedendumd+self.adendumd)

            # clip the involute sides to get rid of overlap

            # new radius, and corresponding x and y

            r = d['r']
            x = r * cos(d['t'])
            y = r * sin(d['t'])

            dx = self.dx(d['t'])
            dy = self.dy(d['t'])

            dx2 = self.dx2(d['t'])
            dy2 = self.dy2(d['t'])

            dx, dy = normalize(dx,dy)

            # tooth width = circular pitch / 2.0 - backlash
            # offset = tooth width / 2.0 + x(td)

            toothWidth = self.cpitch / 2.0
            offset = (toothWidth / 2.0) + involute(rc, td)[0]

            # make sure the involute won't intersect its symmetric partner

            tx = findClosestDown(inverseInvoluteTableX, offset)
            
            if (tt > tx):
                print "clipping tooth "+`n`+": "+`tt`+" "+`tx`
                tt = tx

            if (dy < 0):
                rot = 3.0*pi/2.0+atan(dx/dy)
            else:
                rot = pi/2.0+atan(dx/dy)

            t = 0.0
            m = 0
            
            #f.write("# tooth #"+`n`+"\n")
            #iterate through and collect points
            pts = []
            while m < self.toothSlices:
                ix, iy = involute(rc, t)
                px = -offset+ix
                py = iy

                npx = px * cos(rot) + py * sin(rot)
                npy = px * -sin(rot) + py * cos(rot)
                
                arr = [(x + npx, y + npy,)]
                
                # right side
                px = offset-ix # flip left-hand side
                py = iy

                npx = px * cos(rot) + py * sin(rot)
                npy = px * -sin(rot) + py * cos(rot)
                
                arr.append((x + npx, y + npy,))
                
                pts.append(arr)
                t += tt / (self.toothSlices-1.0)
                m += 1
            
            for i in range(0, len(pts)):
                f.write("  0\nLINE\n")
                f.write("  8\n %s \n" % ('gear')) #layername
                if(i == 0):
                    #for the first one, we don't have much of a line
                    f.write(" 10\n%f\n" % pts[i][0][0])
                    f.write(" 20\n%f\n" % pts[i][0][1])
                else:
                    f.write(" 10\n%f\n" % pts[i-1][0][0])
                    f.write(" 20\n%f\n" % pts[i-1][0][1])
                f.write(" 11\n%f\n" % pts[i][0][0])
                f.write(" 21\n%f\n" % pts[i][0][1])
                
                f.write("  0\nLINE\n")
                f.write("  8\n %s \n" % ('gear'))
                if(i == 0):
                    f.write(" 10\n%f\n" % pts[i][1][0])
                    f.write(" 20\n%f\n" % pts[i][1][1])
                else:
                    f.write(" 10\n%f\n" % pts[i-1][1][0])
                    f.write(" 20\n%f\n" % pts[i-1][1][1])
                f.write(" 11\n%f\n" % pts[i][1][0])
                f.write(" 21\n%f\n" % pts[i][1][1])
            
            #store the first two points to use to connect the teeth
            teeth_ends.append((pts[0][0][0], pts[0][0][1], pts[0][1][0], pts[0][1][1],))
            
            #connect the last two points of the tip of the spur
            f.write("  0\nLINE\n")
            f.write("  8\n %s \n" % ('gear'))
            f.write(" 10\n%f\n" % pts[-1][0][0])
            f.write(" 20\n%f\n" % pts[-1][0][1])
            f.write(" 11\n%f\n" % pts[-1][1][0])
            f.write(" 21\n%f\n" % pts[-1][1][1])
            
            #TODO need to figure out between tooth lines; maybe not, that should be tied together with the above...

            # tooth face (front)
            #f.write("face fft"+`n`+" (")
            #for z in range(self.toothSlices):
            #    f.write(" fpt"+`n`+"_r"+`z`+" ")
            #for z in range(self.toothSlices-1,-1,-1):
            #    f.write(" fpt"+`n`+"_l"+`z`+" ")
            #f.write(") endface\n")

            # (back)
            #f.write("face bft"+`n`+" (")
            #for z in range(self.toothSlices):
            #    f.write(" bpt"+`n`+"_l"+`z`+" ")
            #for z in range(self.toothSlices-1,-1,-1):
            #    f.write(" bpt"+`n`+"_r"+`z`+" ")
            #f.write(") endface\n")

            # tooth sides
            #for z in range(self.toothSlices-1):
            #    f.write("face soft"+`n`+"_l"+`z`+" ( bpt"+`n`+"_l"+`z`+" fpt"+`n`+"_l"+`z`+" fpt"+`n`+"_l"+`(z+1)`+" bpt"+`n`+"_l"+`(z+1)`+") endface\n")
            #    f.write("face soft"+`n`+"_r"+`z`+" ( bpt"+`n`+"_r"+`z+1`+" fpt"+`n`+"_r"+`z+1`+" fpt"+`n`+"_r"+`z`+" bpt"+`n`+"_r"+`z`+") endface\n")

            # tooth side (end)
            #f.write("face soft"+`n`+"_c ( bpt"+`n`+"_l"+`self.toothSlices-1`+" fpt"+`n`+"_l"+`self.toothSlices-1`+" fpt"+`n`+"_r"+`self.toothSlices-1`+" bpt"+`n`+"_r"+`self.toothSlices-1`+") endface\n")

            # in between teeth (gap), on the outside
            #f.write("face sofg"+`n`+" ( fpt"+`(n+1)%self.teethCount`+"_r0 fpt"+`n`+"_l0 bpt"+`n`+"_l0 bpt"+`(n+1)%self.teethCount`+"_r0) endface\n")
            
            # inner point on inner axle
            ri = self.innerradius(d['t'])
            xi, yi = toCartesian(ri, d['t'])

            # front
            #f.write("point fpi"+`n`+" ("+`xi`+" "+`yi`+" "+`self.depth`+" ) endpoint\n")
            # back
            #f.write("point bpi"+`n`+" ("+`xi`+" "+`yi`+" 0.0 ) endpoint\n")
            #store the inner points to create the inner radius
            inner_pts.append((xi,yi,))

            # inner radius faces
            #f.write("face sif"+`n`+" (bpi"+`n`+" fpi"+`n`+" fpi"+`(n+1)%self.teethCount`+" bpi"+`(n+1)%self.teethCount`+") endface\n")

            # gap to inner radius face
            # front
            #f.write("face ffbi"+`n`+" ( fpi"+`(n+1)%self.teethCount`+" fpi"+`n`+" fpt"+`n`+"_l0 fpt"+`(n+1)%self.teethCount`+"_r0) endface\n")
        
            # back
            #f.write("face fbbi"+`n`+" ( bpi"+`(n+1)%self.teethCount`+" bpt"+`(n+1)%self.teethCount`+"_r0 bpt"+`n`+"_l0 bpi"+`n`+" ) endface\n")
            
            # tooth to inner radius face
            # front
            #f.write("face ffbt"+`n`+" ( fpi"+`n`+" fpt"+`n`+"_r0 fpt"+`n`+"_l0) endface\n")
            
            # back
            #f.write("face fbbt"+`n`+" ( bpi"+`n`+" bpt"+`n`+"_l0 bpt"+`n`+"_r0) endface\n")

            print n,
            
            n += 1
        
        #TODO this is the same process repeated twice, put it in a method and test it
        #connect the inner radius
        for i in range(0, len(inner_pts)):
            f.write("  0\nLINE\n")
            f.write("  8\n %s \n" % ('gear')) #layername
            if(i == 0):
                f.write(" 10\n%f\n" % inner_pts[-1][0])
                f.write(" 20\n%f\n" % inner_pts[-1][1])
            else:
                f.write(" 10\n%f\n" % inner_pts[i-1][0])
                f.write(" 20\n%f\n" % inner_pts[i-1][1])
            f.write(" 11\n%f\n" % inner_pts[i][0])
            f.write(" 21\n%f\n" % inner_pts[i][1])
        
        #connect the teeth
        for i in range(0, len(teeth_ends)):
            f.write("  0\nLINE\n")
            f.write("  8\n %s \n" % ('gear')) #layername
            if(i == 0):
                f.write(" 10\n%f\n" % teeth_ends[-1][0])
                f.write(" 20\n%f\n" % teeth_ends[-1][1])
            else:
                f.write(" 10\n%f\n" % teeth_ends[i-1][0])
                f.write(" 20\n%f\n" % teeth_ends[i-1][1])
            f.write(" 11\n%f\n" % teeth_ends[i][2])
            f.write(" 21\n%f\n" % teeth_ends[i][3])
        
        #f.write("# object\n")

        #print "\nassembling polygons into gear object"

        #f.write("object gear (")
        #for n in range(len(self.teethLoc)):
        #    f.write(" fft"+`n`+" ffbi"+`n`+" ffbt"+`n`+" bft"+`n`+" fbbi"+`n`+" fbbt"+`n`+
        #            " sif"+`n`+" soft"+`n`+"_c sofg"+`n`+" ")
        #    for z in range(self.toothSlices-1):
        #        f.write(" soft"+`n`+"_l"+`z`+" soft"+`n`+"_r"+`z`+" ")
        #f.write(") endobject\n")

        #f.write("### END BASE GEOMETRY ###\n")
        #f.write("### TEETH GEOMETRY ###\n")
        #f.write("### "+`len(self.teethLoc)`+" teeth\n")
        
        #f.write("### END TEETH GEOMETRY ###\n")
        #f.write("### END GEAR GEOMETRY ###\n")
        
    def preamble(self,f):
        f.write("  999\nDXF created from gearsgen.py phill baker\n")
        f.write("  0\nSECTION\n")
        f.write("  2\nENTITIES\n")

    def postamble(self,f):
        f.write("  0\nENDSEC\n")
        f.write("  0\nEOF\n")

    def radiusOfCurvature(self, t):
        dx = self.dx(t)
        dy = self.dy(t)
        dx2 = self.dx2(t)
        dy2 = self.dy2(t)
              
        return ((dx*dx+dy*dy)**1.5) / (dx*dy2-dx2*dy)

    #def perimeter(self):
    #    "Needs to be overridden."
    #    return self.oradius*pi

    def calcPoints(self):
        refinement = 1
        lccd = 0.0

        while 1:
            module = self.cpitch / pi

            self.dedendumd = module * 1.25
            self.adendumd = module

            ccd = 0.0
                
            theta = 0.0
            ro = self.outerradius(0) - self.dedendumd
            
            fx = lx = ro * cos(0)
            fy = ly = ro * sin(0)

            actualTeeth = 1

            while theta < 2.0 * pi:
                theta += self.tolerance
                ro = self.outerradius(theta) - self.dedendumd
                x = ro * cos(theta)
                y = ro * sin(theta)
                # add to distance from last point
                # by adding cartesian distance
                ccd += sqrt((x - lx)**2 + (y-ly)**2)
                if (ccd >= self.cpitch and actualTeeth < self.teethCount):
                    # add the new set of points
                    actualTeeth += 1
                    ccd = 0.0

                lx = x
                ly = y

            print "refinement", refinement, ",", "teeth =", actualTeeth, "extra = ", ccd, "c pitch = ", self.cpitch
            #if pastvalues.has_key((actualTeeth, ccd)):
            if self.teethCount == actualTeeth and \
                   (abs(ccd - self.cpitch) < self.tolerance * 10):
                break

#                pastvalues[(actualTeeth, ccd)] = self.cpitch
            refinement += 1
            self.cpitch = (self.cpitch*actualTeeth + ccd) / (self.teethCount + 1)
            lccd = ccd

        print "new circular pitch =",self.cpitch
        print "new perimeter =",self.cpitch*self.teethCount
        
        module = self.cpitch / pi

        self.dedendumd = module * 1.25
        self.adendumd = module

        print "dedendum distance =", self.dedendumd
        print "adendum distance =", self.adendumd
        print "tooth height =", self.dedendumd+self.adendumd

        # compute final refinement and store data

        ccd = 0.0
        theta = 0.0
        ro = self.outerradius(theta) - self.dedendumd
        
        lx = ro * cos(theta)
        ly = ro * sin(theta)

        rc = self.radiusOfCurvature(theta)

        self.teethLoc = [{'x':lx, 'y':ly, 'r':ro, 't':0.0, 'rc':rc}]

        print "computing points for final refinement"

        while len(self.teethLoc) < self.teethCount:
            theta += self.tolerance
            ro = self.outerradius(theta) - self.dedendumd
            x = ro * cos(theta)
            y = ro * sin(theta)
            # add to distance from last point
            # by adding cartesian distance
            ccd = sqrt((x - lx)**2 + (y-ly)**2)
            if (ccd >= self.cpitch):
                # add the new set of points
                rc = self.radiusOfCurvature(theta)
                self.teethLoc.append({'x':x, 'y':y, 'r':ro, 't':theta, 'rc':rc})
                lx = x
                ly = y
                
if (__name__ == "__main__"):
    s = SLFMaker()
    print "### SLFMake self-test"

class SLFMakerSphereTeeth(SLFMaker):
    def __init__(self, teethCount, depth=0.1, tolerance = 0.001):
        SLFMaker.__init__(self)
        
    def amble(self, f):
        """Takes the pointlist and produces the necessary SLIDE code for
        the points and the line they're in"""
        import string

        f.write("### GEAR GEOMETRY ###\n")
        f.write("### BASE GEOMETRY ###\n")
        f.write("### "+`len(self.teethLoc)`+" actual slices\n")

        f.write("# points\n")

        n = 0
        for d in self.teethLoc:
            f.write("sphere t"+`n`+"\n")
#            f.write("  radius "+`self.cpitch/2.0`+"\n") #"+`d['rc']/5.0`+"\n")
            f.write("  radius "+`d['rc']/10.0`+"\n")
            f.write("endsphere\n")
                    
            n += 1

        n = 0
        
        f.write("group gear\n")
        for d in self.teethLoc:
            f.write("instance t"+`n`+"\n")
            f.write("translate ("+`d['x']`+" " + `d['y']` + " 0)\n")
            f.write("endinstance\n")
            n += 1
            
        f.write("endgroup\n")

        f.write("# faces\n")

        f.write("# object\n")

        f.write("### END BASE GEOMETRY ###\n")
        f.write("### TEETH GEOMETRY ###\n")

        f.write("### END TEETH GEOMETRY ###\n")
        f.write("### END GEAR GEOMETRY ###\n")    
