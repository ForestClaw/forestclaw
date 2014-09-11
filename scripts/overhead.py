#!/usr/bin/env python
import sys

if sys.argv.__len__() < 3:
    print "Usage :  Enter time needed for each portion of the code"
    print "    % overhead <advance> <exchange> <regrid>"
    sys.exit();

a = float(sys.argv[1])
e = float(sys.argv[2])
r = float(sys.argv[3])

o = r + e + a

print " "
# print "%40s %6.1f%%" % ("Overhead as percentage of ADVANCE",100.0*o/a)
# print "%40s %6.1f%%" % ("Overhead as percentage of total",100.0*o/(a + o))
print "%40s %6.1f%%" % ("ADVANCE as percentage of total",100.0*a/o)
print "%40s %6.1f%%" % ("EXCHANGE as percentage of total",100.0*e/o)
print "%40s %6.1f%%" % ("REGRID as percentage of total",100.0*r/o)
print " "

print "Grid advances take %5.2f times as long as exchange and regridding" % (a/o)
print " "
