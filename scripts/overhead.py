#!/usr/bin/env python
import sys

if sys.argv.__len__() < 3:
    print "Usage :  Enter time needed for each portion of the code"
    print "    % overhead <advance> <exchange> <regrid>"
    sys.exit();

a = float(sys.argv[1])
e = float(sys.argv[2])
r = float(sys.argv[3])

o = r + e

print " "
print "%40s %6.1f%%" % ("Over head as percentage of ADVANCE",100.0*o/a)
print "%40s %6.1f%%" % ("Over head as percentage of total",100.0*o/(a + o))
print " "

print "Grid advances take %2.1f times as long as exchange and regridding" % (a/o)
print " "
