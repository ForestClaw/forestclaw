
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 201
    nypoints = 201
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    outfile= "flatbathy.topotype2"     
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 101
    nypoints = 101
    xlower = -100.e0
    xupper = 100.e0
    yupper = 100.e0
    ylower = -100.e0
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    """
    flatbathy
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    z = -100*ones(x.shape)
    return z


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    ze = -((x + 50)**2)/100.
    z = 10*exp(ze)

    return where(z > 1e-3,z,0)

if __name__=='__main__':
    maketopo()
    makeqinit()
