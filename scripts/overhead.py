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
t = r + e + a

print " "
print "%40s %6.1f%%" % ("ADVANCE as percentage of total",100.0*a/t)
print "%40s %6.1f%%" % ("EXCHANGE as percentage of total",100.0*e/t)
print "%40s %6.1f%%" % ("REGRID as percentage of total",100.0*r/t)
print " "

if o > a:
    print "Exchange and regridding take %5.2f times as long as grid advances" % (o/a)
else:
    print "Grid advances take %5.2f times as long as exchange and regridding" % (a/o)
print " "
