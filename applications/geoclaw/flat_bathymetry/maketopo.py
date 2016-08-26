
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
    z = ones(x.shape)*-100.e0
    #z = zeros(x.shape)
    return z


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    # ze = -((x+0e0)**2)/10.
    # z = 1.0 + 0.01*exp(ze)
    ze = -((x+0e0)**2 + (y+0e0)**2)/10.
    z = where(ze>-10., 40.e0*exp(ze), 0.)
    return z

if __name__=='__main__':
    maketopo()
    makeqinit()
