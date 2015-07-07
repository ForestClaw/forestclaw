#!/usr/bin/env python
import sys

if sys.argv.__len__() < 4:
    print "Usage :  Enter time needed for each portion of the code"
    print "    % overhead <walltime> <advance> <exchange> <regrid>"
    sys.exit();

if (sys.argv.__len__() == 4):
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

else:
    w = float(sys.argv[1])
    a = float(sys.argv[2])
    e = float(sys.argv[3])
    r = float(sys.argv[4])

    o = w - (a + e + r)

    print " "
    print "%40s %6.1f%%" % ("ADVANCE as percentage of wall time",100.0*a/w)
    print "%40s %6.1f%%" % ("EXCHANGE as percentage of wall time",100.0*e/w)
    print "%40s %6.1f%%" % ("REGRID as percentage of wall time",100.0*r/w)
    print "%40s %6.1f%%" % ("OTHER as percentage of wall time",100.0*o/w)
    print " "
