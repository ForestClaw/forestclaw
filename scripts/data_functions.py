import parallel_tools
import numpy as np
import sys


def fraction_amr(job,mx=None,proc=None,level=None,all=None):

    w = job["walltime"]
    p = job["partition_comm"]
    r = job["regrid"] + job["adapt_comm"]
    i = job["init"]
    g = job["ghostfill"] + job["ghostpatch_comm"]
    c = job["cfl_comm"]
    cost_per_grid = (i+r+g+c+p)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_advance(job,mx=None,proc=None,level=None,all=None):

    a = job["advance"]
    w = job["walltime"]
    cost_per_grid = (a)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_regrid(job,mx=None,proc=None,level=None,all=None):

    p = job["partition_comm"]
    r = job["regrid"] + job["adapt_comm"]
    w = job["walltime"]
    if (r) == 0:
        print "fraction_regrid : Zero regrid time"
        print "mx = %d: proc : %d; level = %d" % (mx,proc,level)
        sys.exit()
    cost_per_grid = (r+p)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_init(job,mx=None,proc=None,level=None,all=None):

    i = job["init"]
    w = job["walltime"]
    cost_per_grid = i/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_ghostcomm(job,mx=None,proc=None,level=None,all=None):

    gc = job["ghostpatch_comm"]
    w = job["walltime"]
    cost_per_grid = (gc)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int

def fraction_ghostfill(job,mx=None,proc=None,level=None,all=None):

    gf = job["ghostfill"]
    w = job["walltime"]
    cost_per_grid = (gf)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_advance(job,mx=None,proc=None,level=None,all=None):

    w = job["walltime"]
    a = job["advance"]
    cost_per_grid = (a)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_cfl(job,mx=None,proc=None,level=None,all=None):

    w = job["walltime"]
    cfl = job["cfl_comm"]
    cost_per_grid = (cfl)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int


def fraction_partition(job,mx=None,proc=None,level=None,all=None):

    w = job["walltime"]
    p = job["partition_comm"]
    cost_per_grid = (p)/w
    fmt_int = False

    v = cost_per_grid

    return v, fmt_int

def fraction_unaccounted(job,mx=None,proc=None,level=None,all=None):

    w = job["walltime"]
    a = job["advance"]
    i = job["init"]
    cost_per_grid = (w-a)/w
    v = cost_per_grid

    fmt_int = False

    return v, fmt_int

def fraction_step1(job,mx=None,proc=None,level=None,all=None):

    g = job["ghostfill"]
    s1 = job["step1"]
    cost_per_grid = (s1)/g
    fmt_int = False

    v = cost_per_grid

    return v, fmt_int

def fraction_step2(job,mx=None,proc=None,level=None,all=None):

    g = job["ghostfill"]
    s2 = job["step2"]
    cost_per_grid = (s2)/g
    fmt_int = False

    v = cost_per_grid

    return v, fmt_int

def fraction_step3(job,mx=None,proc=None,level=None,all=None):

    g = job["ghostfill"]
    s3 = job["step3"]
    cost_per_grid = (s3)/g
    fmt_int = False

    v = cost_per_grid

    return v, fmt_int

def fraction_ghostfill_steps(job,mx=None,proc=None,level=None,all=None):

    g = job["ghostfill"]
    s1 = job["step1"]
    s2 = job["step2"]
    s3 = job["step3"]
    cost_per_grid = (s1+s2+s3)/g
    fmt_int = False

    v = cost_per_grid

    return v, fmt_int

def fraction_ghostfill_unaccounted(job,mx=None,proc=None,level=None,all=None):

    g = job["ghostfill"]
    s1 = job["step1"]
    s2 = job["step2"]
    s3 = job["step3"]
    cost_per_grid = (g-s1-s2-s3)/g
    fmt_int = False

    v = cost_per_grid

    return v, fmt_int



def grids_per_proc(job,mx=None,proc=None,level=None,all=None):
    gpp = job["grids_per_proc"]

    v = gpp
    fmt_int = True

    return v,fmt_int

def grids_per_time(job,mx=None,proc=None,level=None,all=None):
    v = job["advance_steps"]/job["walltime"]
    fmt_int = True

    return v,fmt_int

def grids_per_advance_time(job,mx=None,proc=None,level=None,all=None):
    v = job["advance_steps"]/job["advance"]
    fmt_int = True

    return v,fmt_int


def cost_per_grid(job,mx=None,proc=None,level=None,all=None):
    v = job["walltime"]/job["advance_steps"]

    fmt_int = False

    return v,fmt_int
