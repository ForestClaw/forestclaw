import matplotlib.pyplot as plt
import os
import sys
import numpy as np


def plot_results(val2plot = "walltime"):
    # val2plot = sys.argv[1]

    print val2plot

    data = np.loadtxt("results.out")

    grids_per_edge = 2

    cl = 'r'
    sl = 'o'

    job = {}
    job["jobid"] = data[:,0]
    job["procs"] = data[:,1]
    job["mx"] = data[:,2]
    job["nout"] = data[:,3]
    job["grids"] = data[:,4]
    job["walltime"] = data[:,5]
    job["advance"] = data[:,6]
    job["exchange"] = data[:,7]
    job["regrid"] = data[:,8]
    job["ghostcomm"] = data[:,9]
    job["cfl"] = data[:,10]
    job["extra"] = data[:,11]
    job["time_per_grid"] = data[:,12]

    procs_unique = np.unique(job["procs"])
    njobs = data.shape[0]

    if val2plot in ['walltime','advance','exchange','time_per_grid']:
        y_all = job[val2plot]
    elif val2plot == 'regrid':
        y_all = (job["regrid"] - job["extra"])/grids_advanced
    else:
        y_all = job[val2plot]

    y_avg = []
    y_err = []
    for p in procs_unique:
        y_unique = y_all[np.where(job["procs"] == p)]
        y_avg.append(np.average(y_unique))
        y_err.append(np.std(y_unique))

    plt.errorbar(procs_unique,y_avg,yerr = y_err,fmt='o',hold=True)
    ax = plt.gca();
    ax.set_yscale('log')
    plt.loglog(procs_unique,y_avg,cl+sl+'-',markersize=10,hold=True)
    # plt.loglog(job["procs"],y_all,'b'+sl,markersize=4,hold=True)

    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))
    ax.xaxis.set_major_locator(plt.FixedLocator(procs_unique))
    plt.setp(ax.get_xticklabels(),fontsize=14)
    plt.setp(ax.get_yticklabels(),fontsize=14)
    ax.set_xlabel("Processor count",fontsize=16)
    ax.set_ylabel("%s" % (val2plot.capitalize()),fontsize=16)

    plt.xlim([4**(-0.5), 4**(3.4)])
    plt.ylim([np.min(y_all), np.max(y_all)])
    plt.title(val2plot.capitalize())

    plt.savefig("%s.png" % val2plot)
    plt.show()


if __name__ == "__main__":
    import sys
    if sys.argv.__len__() < 2:
        print "Usage :  Enter one of \"walltime\", \"advance\", \"exchange\"" + \
        "\"regrid\", \"ghostcomm\""
        print "Default : walltime"
        plot_results()
        # sys.exit()
    else:
        plot_results(sys.argv[1])
