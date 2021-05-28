#!/usr/bin/env python3
import os
import re

srcdir = os.path.join(os.getcwd(),'src')
appdir = os.path.join(os.getcwd(),'applications')
testdir = os.path.join(os.getcwd(),'test')

num_errors = 0
for rootdir in [srcdir, appdir, testdir]:
    for subdir, dirs, files in os.walk(rootdir):
        has_makefile = False
        has_cmakefile = False
        named_files = {}
        for file in files:
            has_makefile |= (file=="Makefile.am")
            has_cmakefile |= (file=="CMakeLists.txt")
            if re.match(".*\.apps",file):
                name = file[0:-5]
                if name in named_files:
                    named_files[name] |= {'apps'}
                else:
                    named_files[name] = {'apps'}
            if re.match(".*\.cmake",file):
                name = file[0:-6]
                if name in named_files:
                    named_files[name] |= {'cmake'}
                else:
                    named_files[name] = {'cmake'}
        if has_cmakefile and not has_makefile:
            print(subdir + "\n\t has a CMakeLists.txt file but no Makefile.am file\n")
            num_errors+=1
        if not has_cmakefile and has_makefile:
            print(subdir + "\n\t has a Makefile.am file but no CMakeLists.txt file\n")
            num_errors+=1
        for name in named_files:
            files = named_files[name]
            if 'cmake' in files and 'apps' not in files:
                print(subdir + "\n\t has a " + name + ".cmake file but no " + name +".apps file\n")
                num_errors+=1
            if 'cmake' not in files and 'apps' in files:
                print(subdir + "\n\t has a " + name + ".apps file but no " + name +".cmake file\n")
                num_errors+=1


exit(num_errors)
