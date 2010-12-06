#
#  Generate DXF files for involute gears
#
#  Copyright (C) 2009  Clifford Wolf <clifford@clifford.at>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  ---------------------------------------------------------------------------
#
#  Run 'python gearsgen.py --help' for instructions.
#

import getopt, sys
from math import *

# internal gear
internal = 0

# number of teeth
N = 5

# module
m = 1.0

# clearance
c = 0.0

# addendum offset
off_a = 0.0

# addendum/dedendum multiplier
mul_a = 0.5
mul_d = 0.5

# pitch
p = 0.0

# spacing
s = 0.0

# tip relief
tr = 0.0

# x- and y-offset, rotation
off_x, off_y, off_r = 0.0, 0.0, 0.0

# skip dxf header/footer generation
e_mode = 0

# layer name for generated dxf data
layername = ""

# output filename
outfile = ""

def print_help():
    print("")
    print("Usage: python {0} [ options ]".format(sys.argv[0]))
    print("")
    print("")
    print("Options:")
    print("")
    print(" -i ....... generate outline for internal gear")
    print("            (i.e. the meaning of pitch and clearance is inverted,")
    print("            the sign of the spacing is changed and the tip relief")
    print("            is reset to 0.0)")
    print("")
    print(" -n N ..... number of teeth (default = 5)")
    print("")
    print(" -m U ..... module (default = 1.0)")
    print("            (for gears to mesh, their modules must be equal. the gears")
    print("            working diameter is [number of teeth] times [module])")
    print("")
    print(" -c X ..... clearance (default = 0.0)")
    print("            (move the gear bottom land down by this number of milimeters)")
    print("")
    print(" -a X ..... addendum offset (default = 0.0)")
    print("            (move the gear top land up by this number of milimeters by")
    print("            incrasing the addendum by this value. use with care!)")
    print("")
    print(" -p X ..... pitch (default = 0.0)")
    print("            (move the gear top land up by this number of milimeters by")
    print("            linearly extending the top land. do not confuse with '-a'!)")
    print("")
    print(" -s X ..... spacing (default = 0.0)")
    print("            (move the flanks inside by this number of milimeters.)")
    print("")
    print(" -t U ..... tip relief (default = 0.0)")
    print("            (cut off a triangle with this side length from the gear tooths)")
    print("")
    print(" -D U ..... dedendum multiplier (default = 0.5)")
    print(" -A U ..... addendum multiplier (default = 0.5)")
    print("            (multiply module with this value to calculate dedendum/addendum)")
    print("")
    print(" -e ....... generate only DXF entities without DXF header or DXF footer")
    print("            (usefull for creating larger dxf files from a script)")
    print("")
    print(" -l T ..... name of the DXF layer for generated entities")
    print("            (without this option, no layer information is added to the DXF)")
    print("")
    print(" -x X ..... add this value to all X-coordinates in the generated DXF file")
    print(" -y X ..... add this value to all Y-coordinates in the generated DXF file")
    print(" -r X ..... rotate the generate wheel by this ange in degrees")
    print("            (this can be used together with -e and -l to generate more complex")
    print("            DXF files from scripts. default oriantation is gear at 0/0 with")
    print("            a bottom land facing to the positive x-axis)")
    print("")
    print(" -f T ..... name of the DXF file to generate")
    print("            (without this option, a (tk) window is opend and the generated gear")
    print("            is displayed. note that the gear is scaled to fit the viewport.)")
    print("")
    print("")
    print("Argument types as used in above options reference:")
    print("")
    print(" N = positive integer")
    print(" U = unsigned decimal number")
    print(" X = signed decimal number")
    print(" T = text")
    print("")

try:
    opts, args = getopt.getopt(sys.argv[1:], "in:m:c:a:p:s:t:D:A:el:x:y:r:f:")
except getopt.GetoptError:
    print_help()
    exit(1)

if len(args) > 0:
    print_help()
    exit(1)

for o, a in opts:
    if o == "-i":
        internal = 1
    if o == "-n":
        N = int(a)
    elif o == "-m":
        m = float(a)
    elif o == "-c":
        c = float(a)
    elif o == "-a":
        off_a = float(a)
    elif o == "-p":
        p = float(a)
    elif o == "-s":
        s = float(a)
    elif o == "-t":
        tr = float(a)
    elif o == "-D":
        mul_d = float(a)
    elif o == "-A":
        mul_a = float(a)
    elif o == "-e":
        e_mode = 1
    elif o == "-l":
        layername = a
    elif o == "-x":
        off_x = float(a)
    elif o == "-y":
        off_y = float(a)
    elif o == "-r":
        off_r = float(a)
    elif o == "-f":
        outfile = a

linedata = [];

def line(x1, y1, x2, y2):
    global linedata
    linedata.append(x1)
    linedata.append(y1)
    linedata.append(x2)
    linedata.append(y2)

def involute_intersect_with_r2(r1, r2, phi, phi_step):
    while 1:
        x = r1 * phi
        sx = cos(phi) * r1
        sy = sin(phi) * r1
        px = sx + sin(phi) * x;
        py = sy - cos(phi) * x;
        if sqrt(px**2 + py**2) > r2:
            if phi_step > 1e-10:
                return involute_intersect_with_r2(r1, r2, phi-phi_step, phi_step/2)
            return [ phi, atan2(py, px) ]
        phi = phi + phi_step

def involute(r, phi_start, phi_stop, steps):
    ox = cos(phi_start) * r;
    oy = sin(phi_start) * r;
    for i in range(1, steps+1):
        phi = phi_start + (phi_stop - phi_start)*i/steps
        x = r * (phi - phi_start)
        sx = cos(phi) * r
        sy = sin(phi) * r
        px = sx + sin(phi) * x;
        py = sy - cos(phi) * x;
        line(ox, oy, px, py)
        ox, oy = px, py
    return [ ox, oy ];

if internal:
    p, c = c, p
    s = s * -1.0
    t = 0.0

Dw = float(m*N)
Dr = float(Dw - 2.0*mul_d*m)
Do = float(Dw + 2.0*(mul_a*m + off_a))

# phi_o = angle to the intersection point of the involute and the Do circle
phi_o = involute_intersect_with_r2(Dr/2, Do/2, 1.0, 1.0)[1]

# phi_tr = angle where the involute intersects the Do-2tr circle
phi_tr = involute_intersect_with_r2(Dr/2, Do/2 - tr, 1.0, 1.0)[0]

# ti = angle to the intersection point of the involute and the Dw circle
ti = involute_intersect_with_r2(Dr/2, Dw/2, 1.0, 1.0)[1]

# ts = angle for implementing spacing at Dw
ts = pi*s / Dw

# t = angle span of bottom land
t = pi/N - ti*2 + ts*2

for i in range(N):
    phi1 = 2*pi*i/N - t/2 + off_r*(pi/180)
    phi2 = phi1 + t
    phi3 = 2*pi*(i+1)/N - t/2 + off_r*(pi/180)
    if c != 0:
        line(cos(phi1)*(Dr/2), sin(phi1)*(Dr/2), cos(phi1)*(Dr/2-c), sin(phi1)*(Dr/2-c))
    line(cos(phi1)*(Dr/2-c), sin(phi1)*(Dr/2-c), cos(phi2)*(Dr/2-c), sin(phi2)*(Dr/2-c))
    if c != 0:
        line(cos(phi2)*(Dr/2), sin(phi2)*(Dr/2), cos(phi2)*(Dr/2-c), sin(phi2)*(Dr/2-c))
    ax, ay = involute(Dr/2, phi2, phi2 + phi_tr, 10)
    bx, by = involute(Dr/2, phi3, phi3 - phi_tr, 10)
    if t > 0:
        axt, ayt = ax, ay
        bxt, byt = bx, by
        phi_xo = phi_o + (tr/2)/(Do/2)
        ax = cos(phi2 + phi_xo) * Do/2
        ay = sin(phi2 + phi_xo) * Do/2
        bx = cos(phi3 - phi_xo) * Do/2
        by = sin(phi3 - phi_xo) * Do/2
        line(axt, ayt, ax, ay)
        line(bxt, byt, bx, by)
    if p != 0:
        aa = atan2(ay, ax)
        ar = sqrt(ay**2 + ax**2) + p
        ba = atan2(by, bx)
        br = sqrt(by**2 + bx**2) + p
        line(ax, ay, cos(aa)*ar, sin(aa)*ar)
        line(cos(aa)*ar, sin(aa)*ar, cos(ba)*br, sin(ba)*br)
        line(bx, by, cos(ba)*br, sin(ba)*br)
    else:
        line(ax, ay, bx, by)

if outfile != "":
    if outfile != "-":
        f = open(outfile, "w")
    else:
        f = sys.stdout

    if not e_mode:
        f.write("  0\nSECTION\n")
        f.write("  2\nENTITIES\n")

    for i in range(len(linedata)//4):
        x1 = linedata[i*4 + 0]
        y1 = linedata[i*4 + 1]
        x2 = linedata[i*4 + 2]
        y2 = linedata[i*4 + 3]
        f.write("  0\nLINE\n")
        if layername != "":
            f.write("  8\n{0}\n".format(layername))
        f.write(" 10\n{0}\n".format(x1 + off_x))
        f.write(" 20\n{0}\n".format(y1 + off_y))
        f.write(" 11\n{0}\n".format(x2 + off_x))
        f.write(" 21\n{0}\n".format(y2 + off_y))

    if not e_mode:
        f.write("  0\nENDSEC\n")
        f.write("  0\nEOF\n")

    if outfile != "-":
        f.close()

else:
    import Tkinter
    tk = Tkinter.Tk()
    tk.wm_geometry("400x400")

    canvas = Tkinter.Canvas(tk, width=400, height=400)
    canvas.pack()

    max_xy = 0

    for i in range(len(linedata)//4):
        max_xy = max(max_xy, abs(linedata[i*4 + 0]))

    scale = 180 / max_xy

    for i in range(len(linedata)//4):
        x1 = linedata[i*4 + 0]
        y1 = linedata[i*4 + 1]
        x2 = linedata[i*4 + 2]
        y2 = linedata[i*4 + 3]
        canvas.create_line(200 + x1*scale, 200 - y1*scale, 200 + x2*scale, 200 - y2*scale)

    x = scale * Dw/2
    canvas.create_oval(200 - x, 200 - x, 200 + x, 200 + x, outline="green")

    tk.mainloop()

exit(0)

