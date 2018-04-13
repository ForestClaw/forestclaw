import os
import sys
import numpy as np
from string import Template
import re
from pdb import set_trace


def compile_results(results_dir=None,results_file='results.out',
                    execname=None,duplicate=False):

    if results_dir == None:
        results_dir = os.getcwd()

    # Stuff to grep for
    stats_list_float = [ 'INIT',
                         'ADVANCE$',
                         'GHOSTFILL$',
                         'REGRID$',
                         'GHOSTPATCH_COMM',
                         'ADAPT_COMM',
                         'PARTITION_COMM',
                         'CFL_COMM',
                         'WALLTIME',
                         'GHOSTFILL_COPY',
                         'LOCAL',
                         'COMM',
                         'PARTITION',
                         'GHOSTPATCH_BUILD',
                         'GHOSTFILL_STEP1',
                         'GHOSTFILL_STEP2',
                         'GHOSTFILL_STEP3',
                         'GHOSTFILL_COPY',
                         'GHOSTFILL_INTERP',
                         'GHOSTFILL_AVERAGE',
                         'LOCAL_BOUNDARY_RATIO',
                         'REMOTE_BOUNDARY_RATIO']

    stats_list_int = ['ADVANCE_STEPS_COUNTER$',
                      'GRIDS_PER_PROC$',
                      'LOCAL_BOUNDARY',
                      'REMOTE_BOUNDARY',
                      'INTERIOR']

    # ------------------------------
    # Options : Variable field width
    # ------------------------------
    option_list = {'jobid'  :('jobid', '{jobid[0]:>6s}',   '{:>6d}'),
                   'p'      :('p',     '{p[0]:>8s}',       '{:>8d}'),
                   'mx'     :('mx',    '{mx[0]:>6s}',      '{:>6d}'),
                   'mi'     :('mi',    '{mi[0]:>4s}',      '{:>4d}'),
                   'mj'     :('mj',    '{mj[0]:>4s}',      '{:>4d}'),
                   'min'    :('min',   '{min[0]:>6s}',     '{:>6d}'),
                   'max'    :('max',   '{max[0]:>6s}',     '{:>6d}'),
                   'nout'   :('nout',  '{nout[0]:>8s}',    '{:>8d}'),
                   'tfinal' :('tfinal','{tfinal[0]:>12s}', '{:>12.2e}')}


    option_list_header = "{jobid[1]:>6s}{p[1]:>8s}{mx[1]:>6s}{mi[1]:>4s}{mj[1]:>4s}" + \
                         "{min[1]:>6s}{max[1]:>6s}{nout[1]:>8s}{tfinal[1]:>12s}"
    option_fmt_header = option_list_header.format(**option_list)  # Process names

    option_str_header = option_fmt_header.format(**option_list)
    option_width = 60     # add up all fields widths, i.e. 6+8+6+4+... = 60

    # ---------------------
    # Ints : Field width 12
    # ---------------------
    int_list = ('adv_steps','grids_proc','local_bdry','remote_bdry','interior')
    int_width = 12*len(int_list)
    int_str = ("{:>12s}"*len(int_list)).format(*int_list)

    # -----------------------
    # Floats : Field width 12
    # -----------------------
    float_list = ('init','advance','ghostfill','regrid','patch_comm',
                  'adapt','partition','cfl','walltime',
                  'gf_copy','local','comm','partition','gbuild','step1',
                  'step2','step3', 'copy',
                  'interp','average','l_ratio','r_ratio')
    float_width = 12*len(float_list)
    float_str = ("{:>12s}"*len(float_list)).format(*float_list)


    # Compile data to this file
    line_len = option_width + int_width + float_width
    header_str = option_str_header + int_str + float_str

    # Write out header to file
    resultsfile = open(results_file,'w')
    # resultsfile.write("# " + "-"*line_len)
    # resultsfile.write("\n")

    resultsfile.write(header_str)

    # resultsfile.write("# " + "-"*line_len)
    resultsfile.write("\n")

    # Look for files that look like "torus_00004.o4567".  But if execname is
    # not supplied, be happy with "*_00004.o4567"
    if not execname == None:
        pattern = re.compile("%s_[0-9]{5}\.o[0-9]*" % execname)
    else:
        pattern = re.compile(".*_[0-9]{5}\.o[0-9]*")

    output_files = os.listdir(results_dir)
    for f in output_files:
        if re.match(pattern,f):
            # Extract proc count and jobid from file name, e.g. torus_00016.o45678
            s = f.partition('_')[2].partition('.o')

            pcount = int(s[0])
            jobid = int(s[2])

            run_file = open(f,'r')
            lines = run_file.readlines()

            found_bad_file = False

            for i,l in enumerate(lines):
                if re.search("exceeded",l):
                    found_bad_file = True
                    break

            #if found_bad_file:
            #    print "Skipping file %s (time exceeded)" % f
            #    continue

            # jobid
            resultsfile.write("  ")   # Make up for "# " added as comments above
            s = "{jobid[2]:s}".format(**option_list)
            resultsfile.write(s.format(jobid))

            # proc count
            resultsfile.write("%8d" % pcount)

            # Get values of options

            # mx
            for i,l in enumerate(lines):
                found_mx = False
                if re.search("mx",l):
                    l1 = lines[i].split()
                    mx = int(l1[2])
                    s = "{mx[2]:s}".format(**option_list)
                    resultsfile.write(s.format(mx))
                    found_mx = True
                    break
                if not found_mx:
                    if re.search("nxmax",l):
                        l1 = lines[i].split()
                        mx = int(l1[2])
                        s = "{mx[2]:s}".format(**option_list)  # for ash3d model
                        resultsfile.write(s.format(mx))
                        found_mx = True
                        break


            # mi (look for something distinctive;  count down from there)
            mi = 1
            mj = 1
            for i,l in enumerate(lines):
                if re.search("manifold",l):  # 'mi' not distinctive enough
                    if not re.search("mi",lines[i+1]):
                        print("Brick count mi not found in {:s}; setting mi=1".format(f))
                        mi = 1
                    else:
                        mi_line = lines[i+1].split()
                        mi = int(mi_line[2])
                    if not re.search("mj",lines[i+2]):
                        print("Brick count mj not found in {:s}; setting mj=1".format(f))
                        mj = 1
                    else:
                        mj_line = lines[i+2].split()
                        mj = int(mj_line[2])

                    s = "{mi[2]:s}{mj[2]:s}".format(**option_list)
                    resultsfile.write(s.format(mi,mj))
                    break
                    # print "mi = %d; mj = %d" %(mi,mj)

            # minlevel
            for i,l in enumerate(lines):
                if re.search("minlevel",l):
                    l1 = lines[i].split()
                    if duplicate:
                        minlevel = int(l1[2]) + np.log(float(mi))/np.log(2.0)
                    else:
                        minlevel = int(l1[2])

                    s = "{min[2]:s}".format(**option_list)
                    resultsfile.write(s.format(int(minlevel)))
                    break

            # maxlevel
            for i,l in enumerate(lines):
                if re.search("[^-]maxlevel",l):  # don't match maxlevel-uses-outstyle
                    l1 = lines[i].split()
                    if duplicate:
                        maxlevel = int(l1[2]) + np.log(float(mi))/np.log(2.0)
                    else:
                        maxlevel = int(l1[2])

                    s = "{max[2]:s}".format(**option_list)
                    resultsfile.write(s.format(int(maxlevel)))
                    break

            # nout
            for i,l in enumerate(lines):
                if re.search("nout",l):
                    l1 = lines[i].split()
                    nout = int(l1[2])
                    s = "{nout[2]:s}".format(**option_list)
                    resultsfile.write(s.format(int(l1[2])))
                    break

            # tfinal
            for i,l in enumerate(lines):
                if re.search("tfinal",l):
                    l1 = lines[i].split()
                    tfinal = float(l1[2])
                    s = "{tfinal[2]:s}".format(**option_list)
                    resultsfile.write(s.format(float(l1[2])))
                    break


            # Get counter values
            for k,w in enumerate(stats_list_int):
                found = False
                for i,l in enumerate(lines):
                    if re.search(w,l):
                        a = i
                        found = True
                        break

                if found:
                    l2 = lines[a+2].split()
                    resultsfile.write("%12d" % np.round(float(l2[5])).astype(int))
                else:
                    print("{:s} in file {:s} not found".format(w,f))
                    resultsfile.write("%12f" % (np.nan))

            # Get timer values
            for k,w in enumerate(stats_list_float):
                found = False
                for i,l in enumerate(lines):
                    if re.search(w,l):
                        a = i
                        found = True
                        break

                if found:
                    l2 = lines[a+2].split()
                    resultsfile.write("%12.4e" % float(l2[5]))
                else:
                    print("{:s} in file {:s} not found".format(w,f))
                    resultsfile.write("%12f" % (np.nan))

            resultsfile.write('\n')

            run_file.close()

    resultsfile.close()
