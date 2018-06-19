import os
import sys
import difflib
import subprocess


geo_files = os.listdir(os.getcwd())

p = '/Users/calhoun/clawpack5.0/clawpack/geoclaw/src/2d/shallow'

diff = []
nodiff = []
notfound = []
for f_fclaw in geo_files:
    f_geoclaw = os.path.join(p,f_fclaw)
    if not os.path.isfile(f_geoclaw):
        notfound += [f_fclaw]
        continue

    f1 = open(f_fclaw,'r')    
    f2 = open(f_geoclaw,'r')

    # Do a standard diff on two files
    arg_list = ['diff', f_fclaw, f_geoclaw]
    pp = subprocess.run(arg_list,stdout=subprocess.PIPE)    


    if pp.returncode > 0:
        diff += [f_fclaw]
    else:
        nodiff += [f_fclaw]

print("\n")
print("Files with differences")
print("----------------------")
for f in diff:
    print("{:s}".format(f))

print("\n")
print("Files with no differences")
print("-------------------------")
for f in nodiff:
    print("{:s}".format(f))

print("\n")
print("Files not found")
print("---------------")
for f in notfound:
    print("{:s}".format(f))



