import regions_utility

dim = 2
regiondata = regions_utility.RegionData(dim)

# Format : [minlevel,maxlevel,t1,t2,xlower,xupper,ylower,yupper]
#  minlevel,maxlevel   : levels over which region is active
#  xlower, xupper,ylower,yupper   : Location of the region
#  t0,t1 : (t0,t1) interval over which region is active.

regiondata.regions.append([0, 0, 0., 1.e10,    -1, 1, -1,1.])

regiondata.regions.append([0, 2, 0., 1  ,       0, 1,    0.,1.])

# Only refine in upper right corner
# regiondata.regions.append([0, 2, 1., 1e10  ,       -1,0.,    -1,0])

regiondata.write(data_source='write_regions.py')


