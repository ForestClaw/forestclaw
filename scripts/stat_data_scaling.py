import numpy as np

import os
import sys
import re
import ConfigParser

config = ConfigParser.SafeConfigParser(allow_no_value=True)
config.read('create_run.ini')

execname  = config.get('Problem','execname').partition('#')[0].strip()
njobs     = int(config.get('Run','njobs').partition('#')[0].strip())
proc0     = int(config.get('Run','proc').partition('#')[0].strip())
mode      = config.get('Run','mode').partition('#')[0].strip()

R = np.array(range(0,njobs))
procs = proc0*4**R

# ----------------------------------------------------------------------------------------
stats_list = [ 'WALLTIME', 'ADVANCE', 'Statistics for EXCHANGE','Statistics for REGRID$',
               'GHOSTCOMM$','Statistics for CFL','Statistics for EXTRA4$','TIME_PER_GRID']
# ----------------------------------------------------------------------------------------

output_files = os.listdir(os.getcwd())

resultsfile = open('results.out','w')
resultsfile.write("# " + "-"*138)
resultsfile.write("\n")
fstr = "# %8s%8s%6s%8s%12s" + "%12s"*8 + "\n"
resultsfile.write(fstr % ('jobid','p','mx','nout','grids/proc','Wall','Advance',
                          'Exchange','Regrid','Comm','cfl','extra','time/grid'))

resultsfile.write("# " + "-"*138)
resultsfile.write("\n")
for j,p in enumerate(procs):
    fname = '%s_%05d' % (execname,p)

    for f in output_files:
        if (f.startswith(fname)):

            # Extract proc count and jobid from output file name
            # E.g. torus_01.o4578
            s = f.split('_')[1].split('.o')
            pcount = procs[j]
            jobid = int(s[1])

            run_file = open(f,'r')
            lines = run_file.readlines()

            # jobid
            resultsfile.write("  ")   # Make up for "# " added as comments above
            resultsfile.write("%8d" % jobid)

            # proc count
            resultsfile.write("%8d" % pcount)

            # mx
            for i,l in enumerate(lines):
                if re.search("mx",l):
                    l1 = lines[i].split()
                    mx = int(l1[2])
                    resultsfile.write("%6s" % l1[2])
                    break

            # nout
            for i,l in enumerate(lines):
                if re.search("nout",l):
                    l1 = lines[i].split()
                    nout = int(l1[2])
                    resultsfile.write("%8s" % l1[2])
                    break

            # minlevel
            for i,l in enumerate(lines):
                if re.search("minlevel",l):
                    l1 = lines[i].split()
                    minlevel = int(l1[2])
                    break

            # maxlevel
            for i,l in enumerate(lines):
                if re.search("maxlevel",l):
                    l1 = lines[i].split()
                    maxlevel = int(l1[2])
                    break

            # Grids per proc (average)
            grids_advanced_uniform = (2**maxlevel)**2/pcount*nout*2**(maxlevel-minlevel)
            if mode == 'adapt':
                for i,l in enumerate(lines):
                    if re.search('TIME_PER_GRID',l):
                        a = i
                        break
                l2 = lines[a+2].split()
                tpg = float(l2[5]);

                for i,l in enumerate(lines):
                    if re.search('WALLTIME',l):
                        a = i
                        break

                l2 = lines[a+2].split()
                wt = float(l2[5])

                grids_advanced_adaptive = float(wt/tpg)
                resultsfile.write("%12.2e" % float(grids_advanced_adaptive/grids_advanced_uniform))

            else:
                resultsfile.write("%12d" % int(grids_advanced_uniform/nout))


            # Get everything else in list
            for k,w in enumerate(stats_list):
                for i,l in enumerate(lines):
                    if re.search(w,l):
                        a = i
                        break

                l2 = lines[a+2].split()
                resultsfile.write("%12.4e" % float(l2[5]))

            resultsfile.write('\n')


            run_file.close()

resultsfile.close()
