#!/usr/bin/python

import os;
import subprocess;
import sys
from numpy import *


mthlim_vec = array([0, 4]);
mthlim_str = array(["un","sp"]);

execfile = 'compute_diag'


for m in mthlim_vec:
    diagfile_out = "diagfile.out" + ".mthlim" + str(m);
    f_diag = open(diagfile_out,'w')
    f_diag.write('')
    f_diag.close()

    if m == 0:
        angle_vec = array([1.5, 0.75, 0.1875]);
        mx_vec = array([120, 240, 960]);
    else:
        angle_vec = array([1.5, 0.75, 0.28125]);
        mx_vec = array([120, 240, 640]);

    ncnt = 1;
    for i in range(0,3):
        a = angle_vec[i];
        mx = mx_vec[i];

        s1 = "ccb-" + str(a) + "-mthlim" + str(m);
        os.chdir(s1);

        print "In %s" % (s1)
        os.system('rm -f compute_mixing');
        os.system('rm -f diag.dat');
        os.system('ln -s ../compute_mixing');
        os.system('ln -s ../diag.dat');

        os.system("compute_mixing");
        os.chdir("..");

        diagfile_in = s1 + "/diagfile.in";
        try:
            f = open(diagfile_in)
            line = f.readline()
            print "%5d (my = %3d) %80s" % (ncnt,mx/2,line.strip())
            f.close()

            f_diag = open(diagfile_out,'a')
            f_diag.write(line)
            f_diag.close()

        except:
            print "File %s not found" % (errfile_in)

        ncnt = ncnt + 1

    print " "
