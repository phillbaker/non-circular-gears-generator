#!/usr/bin/python2.2

from math import sin, cos, pi, sqrt, atan
from slfmaker import *

class OffsetPairSLFMaker:
    def __init__(self, teethCount, ts=5, depth=0.1, tolerance = 0.001, secondOffset=1.2):
        self.depth = depth
        self.secondOffset = secondOffset
        self.teethLoc = []
        self.gapLoc = []
        self.teethCount = teethCount
        self.cpitch = self.perimeter() / self.teethCount
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
        self.doShape(f, self.gapLoc, "a", "sBlu")
        self.doShape(f, self.teethLoc, "b", "sRed")

    def doShape(self, f, teethLoc, gearLetter, color):
        """Takes the pointlist and produces the necessary SLIDE code for
        the points and the line they're in"""
        import string

        f.write("### GEAR "+gearLetter+" GEOMETRY ###\n")
        f.write("# points\n")

        n = 0
        inc = 0.01

        print "teeth completed: ",
        
        for d in teethLoc:
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
            
            if (tt > tx): tt = tx

            if (dy < 0):
                rot = 3.0*pi/2.0+atan(dx/dy)
            else:
                rot = pi/2.0+atan(dx/dy)

            t = 0.0
            m = 0
            f.write("# tooth #"+`n`+"\n")
            while m < self.toothSlices:
                ix, iy = involute(rc, t)
                px = -offset+ix
                py = iy

                npx = px * cos(rot) + py * sin(rot)
                npy = px * -sin(rot) + py * cos(rot)
                
                # left side (front)
                f.write("point "+gearLetter+"fpt"+`n`+"_l"+`m`+" ("+ \
                        `x+npx`+" " + \
                        `y+npy`+" " + \
                        `self.depth`+") endpoint\n")

                # (back)

                f.write("point "+gearLetter+"bpt"+`n`+"_l"+`m`+" ("+ \
                        `x+npx`+" " + \
                        `y+npy`+" " + \
                        "0.0) endpoint\n")
                # right side
                
                px = offset-ix # flip left-hand side
                py = iy

                npx = px * cos(rot) + py * sin(rot)
                npy = px * -sin(rot) + py * cos(rot)

                # (front)
                
                f.write("point "+gearLetter+"fpt"+`n`+"_r"+`m`+" ("+ \
                        `x+npx`+" " + \
                        `y+npy`+" " + \
                        `self.depth`+") endpoint\n")

                # (back)

                f.write("point "+gearLetter+"bpt"+`n`+"_r"+`m`+" ("+ \
                        `x+npx`+" " + \
                        `y+npy`+" " + \
                        "0.0) endpoint\n")
                t += tt / (self.toothSlices-1.0)
                m += 1

            # tooth face (front)
                
            f.write("face "+gearLetter+"fft"+`n`+" (")
            for z in range(self.toothSlices):
                f.write(" "+gearLetter+"fpt"+`n`+"_r"+`z`+" ")
            for z in range(self.toothSlices-1,-1,-1):
                f.write(" "+gearLetter+"fpt"+`n`+"_l"+`z`+" ")
            f.write(") endface\n")

            # (back)

            f.write("face "+gearLetter+"bft"+`n`+" (")
            for z in range(self.toothSlices):
                f.write(" "+gearLetter+"bpt"+`n`+"_l"+`z`+" ")
            for z in range(self.toothSlices-1,-1,-1):
                f.write(" "+gearLetter+"bpt"+`n`+"_r"+`z`+" ")
            f.write(") endface\n")

            # tooth sides

            for z in range(self.toothSlices-1):
                f.write("face "+gearLetter+"soft"+`n`+"_l"+`z`+" ( "+gearLetter+"bpt"+`n`+"_l"+`z`+" "+gearLetter+"fpt"+`n`+"_l"+`z`+" "+gearLetter+"fpt"+`n`+"_l"+`(z+1)`+" "+gearLetter+"bpt"+`n`+"_l"+`(z+1)`+") endface\n")
                f.write("face "+gearLetter+"soft"+`n`+"_r"+`z`+" ( "+gearLetter+"bpt"+`n`+"_r"+`z+1`+" "+gearLetter+"fpt"+`n`+"_r"+`z+1`+" "+gearLetter+"fpt"+`n`+"_r"+`z`+" "+gearLetter+"bpt"+`n`+"_r"+`z`+") endface\n")

            # tooth side (end)

            f.write("face "+gearLetter+"soft"+`n`+"_c ( "+gearLetter+"bpt"+`n`+"_l"+`self.toothSlices-1`+" "+gearLetter+"fpt"+`n`+"_l"+`self.toothSlices-1`+" "+gearLetter+"fpt"+`n`+"_r"+`self.toothSlices-1`+" "+gearLetter+"bpt"+`n`+"_r"+`self.toothSlices-1`+") endface\n")

            # in between teeth (gap), on the outside

            f.write("face "+gearLetter+"sofg"+`n`+" ( "+gearLetter+"fpt"+`(n+1)%self.teethCount`+"_r0 "+gearLetter+"fpt"+`n`+"_l0 "+gearLetter+"bpt"+`n`+"_l0 "+gearLetter+"bpt"+`(n+1)%self.teethCount`+"_r0) endface\n")
            
            # inner point on inner axle

            ri = self.innerradius(d['t'])
            xi, yi = toCartesian(ri, d['t'])

            # front
            f.write("point "+gearLetter+"fpi"+`n`+" ("+`xi`+" "+`yi`+" "+`self.depth`+" ) endpoint\n")
            # back
            f.write("point "+gearLetter+"bpi"+`n`+" ("+`xi`+" "+`yi`+" 0.0 ) endpoint\n")

            # inner radius faces

            f.write("face "+gearLetter+"sif"+`n`+" ("+gearLetter+"bpi"+`n`+" "+gearLetter+"fpi"+`n`+" "+gearLetter+"fpi"+`(n+1)%self.teethCount`+" "+gearLetter+"bpi"+`(n+1)%self.teethCount`+") endface\n")

            # gap to inner radius face

            # front
            f.write("face "+gearLetter+"ffbi"+`n`+" ( "+gearLetter+"fpi"+`(n+1)%self.teethCount`+" "+gearLetter+"fpi"+`n`+" "+gearLetter+"fpt"+`n`+"_l0 "+gearLetter+"fpt"+`(n+1)%self.teethCount`+"_r0) endface\n")
        
            # back
            f.write("face "+gearLetter+"fbbi"+`n`+" ( "+gearLetter+"bpi"+`(n+1)%self.teethCount`+" "+gearLetter+"bpt"+`(n+1)%self.teethCount`+"_r0 "+gearLetter+"bpt"+`n`+"_l0 "+gearLetter+"bpi"+`n`+" ) endface\n")
            # tooth to inner radius face

            # front
            f.write("face "+gearLetter+"ffbt"+`n`+" ( "+gearLetter+"fpi"+`n`+" "+gearLetter+"fpt"+`n`+"_r0 "+gearLetter+"fpt"+`n`+"_l0) endface\n")
            
            # back
            f.write("face "+gearLetter+"fbbt"+`n`+" ( "+gearLetter+"bpi"+`n`+" "+gearLetter+"bpt"+`n`+"_l0 "+gearLetter+"bpt"+`n`+"_r0) endface\n")

            print n,
            
            n += 1

        f.write("# object\n")

        print "\nassembling polygons into gear"+gearLetter+" object"

        f.write("object gear"+gearLetter+" (")
        for n in range(len(teethLoc)):
            f.write(" "+gearLetter+"fft"+`n`+" "+gearLetter+"ffbi"+`n`+" "+gearLetter+"ffbt"+`n`+" "+gearLetter+"bft"+`n`+" "+gearLetter+"fbbi"+`n`+" "+gearLetter+"fbbt"+`n`+
                    " "+gearLetter+"sif"+`n`+" "+gearLetter+"soft"+`n`+"_c "+gearLetter+"sofg"+`n`+" ")
            for z in range(self.toothSlices-1):
                f.write(" "+gearLetter+"soft"+`n`+"_l"+`z`+" "+gearLetter+"soft"+`n`+"_r"+`z`+" ")
        f.write(")\n")
        f.write("  surface "+color+"\n")
        f.write("endobject\n")

        f.write("### END GEAR "+gearLetter+" GEOMETRY ###\n")
        
    def preamble(self,f):
        f.write("""
######################
# ncgear.slf, generated by SLFMaker.py
# Jeff Schoner
######################

################## INITIALIZATIONS #########################

###  Get some generic capabilities  ###
tclinit { 
  package require slideui
  set rot 0
}

tclupdate {
#  set rot [expr $rot+1]
}


###  Display window  ###
tclinit { 
  toplevel .slfWindow.gRoot
  CreateSLIDEObjectObject gRoot
  set widget [CreateSLIDEGroupUI .slfWindow.gRoot gRoot]
  pack $widget
}
 

### SURFACES ###

surface sRed  color (1 0.3 0.3) endsurface
surface sGrn  color (0.2 1 0.2) endsurface
surface sBlu  color (0.3 0.3 1) endsurface
surface sYel  color (1 0.8 0) endsurface

""")
    def postamble(self,f):
        width = self.width() / 2.0
        f.write("""

### Scene Assembly ###

group gScene
instance geara
    translate ("""+`self.centerOffset-self.secondOffset`+""" 0 0)
endinstance
instance gearb
    translate ("""+`self.centerOffset+self.secondOffset`+""" 0 0)
endinstance
#instance sgear
#endinstance
endgroup

#################### 
# CAMERA
#################### 

camera cam
  projection SLF_PARALLEL
  frustum (-"""+`width`+" -"+`width`+" -100) ("+`width`+" "+`width`+""" 100)
endcamera

group gCamera
  instance cam
    id instCam
    translate ( 0.0 0.0 1 )
  endinstance
endgroup

#################### 
# LIGHT
#################### 

light lite
  type SLF_DIRECTIONAL
endlight

group gLight
  instance lite
    id instLite

    lookat
      eye ( 1.0 1.0 1.0 )
      target ( 0.0 0.0 0.0 )
      up ( 0.0 1.0 0.0 )
    endlookat

  endinstance
endgroup

light lite2
  type SLF_AMBIENT
  color (0.5 0.5 0.5)
endlight

group gLight2
  instance lite2
    id instLite2
  endinstance
endgroup

####################
# RENDER
####################

window Window
  background (0.3 0.6 0.9)
endwindow

viewport vp Window
  origin ( 0.0 0.0 )
  size ( 1.0 1.0 )
endviewport

render vp gCamera.instCam.cam gScene
  light gLight.instLite.lite
  light gLight2.instLite2.lite2
endrender
##########################################
""")

    def radiusOfCurvature(self, t):
        dx = self.dx(t)
        dy = self.dy(t)
        dx2 = self.dx2(t)
        dy2 = self.dy2(t)
              
        return ((dx*dx+dy*dy)**1.5) / (dx*dy2-dx2*dy)

    def calcPoints(self):
        refinement = 1
        halfteethCount = self.teethCount * 2
        halfpitch = self.cpitch / 2.0

        while 1:
            module = halfpitch * 2.0 / pi

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
                if (ccd >= halfpitch and actualTeeth < halfteethCount):
                    # add the new set of points
                    actualTeeth += 1
                    ccd = 0.0

                lx = x
                ly = y

            print "refinement", refinement, ",", "half-teeth =", actualTeeth, "extra = ", ccd, "c pitch = ", halfpitch * 2.0
            #if pastvalues.has_key((actualTeeth, ccd)):
            if halfteethCount == actualTeeth and \
                   (abs(ccd - halfpitch) < self.tolerance * 5):
                break

#                pastvalues[(actualTeeth, ccd)] = self.cpitch
            refinement += 1
            halfpitch = (halfpitch*actualTeeth + ccd) / (halfteethCount + 1)


        self.cpitch = halfpitch * 2.0

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
        self.gapLoc = []
        gapFlag = 1

        print "computing points for final refinement"

        actualTeeth = 1

        while actualTeeth < halfteethCount:
            theta += self.tolerance
            ro = self.outerradius(theta) - self.dedendumd
            x = ro * cos(theta)
            y = ro * sin(theta)
            # add to distance from last point
            # by adding cartesian distance
            ccd = sqrt((x - lx)**2 + (y-ly)**2)
            if (ccd >= halfpitch):
                # add the new set of points
                rc = self.radiusOfCurvature(theta)
                if gapFlag:
                    self.gapLoc.append({'x':x, 'y':y, 'r':ro, 't':theta, 'rc':rc})
                else:
                    self.teethLoc.append({'x':x, 'y':y, 'r':ro, 't':theta, 'rc':rc})

                gapFlag = not gapFlag
                lx = x
                ly = y
                actualTeeth += 1

if (__name__ == "__main__"):
    s = SLFMaker()
    print "### OffsetPairSLFMake self-test"

