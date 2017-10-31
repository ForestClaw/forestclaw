import numpy as np

def read_topo_data(topofile):

    f = open(topofile,'r')

    l = f.readline()
    ncols = np.fromstring(l.split()[0].strip(),sep=' ')                    

    l = f.readline()
    nrows = np.fromstring(l.split()[0].strip(),sep=' ')

    l = f.readline()
    xllcorner = np.fromstring(l.split()[0].strip(),sep=' ')

    l = f.readline()
    yllcorner = np.fromstring(l.split()[0].strip(),sep=' ')

    l = f.readline()
    cellsize = np.fromstring(l.split()[0].strip(),sep=' ')

    return ncols[0],nrows[0],xllcorner[0],yllcorner[0],cellsize[0]

def compute_distances(lon,lat):

    # lon = [ax bx bx ax ax];
    # lat = [ay ay by by ay];

    lon = np.array([lon[0],lon[1],lon[1],lon[0],lon[0]])
    lat = np.array([lat[0],lat[0],lat[1],lat[1],lat[0]])

    labels = ['top', 'right','bottom','left']

    d = np.zeros(4)
    for i in range(0,4):
        lat1 = np.deg2rad(lat[i])
        lat2 = np.deg2rad(lat[i+1])
        lon1 = np.deg2rad(lon[i])
        lon2 = np.deg2rad(lon[i+1])
    
        dlon = lon2-lon1
        dlat = lat1-lat2
    
        R = 6367.5e3
        a = (np.sin(dlat/2.0))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2.0))**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        d[i] = R * c

        print("%20s %12.4f" % (labels[i],d[i]))

    print("")
    avg_width = (d[0] + d[2])/2.0
    avg_height = (d[1] + d[3])/2.0
    print("%20s %12.4f" % ("Average width",avg_width))
    print("%20s %12.4f" % ("Average height",avg_height))

    print("")
    print("%20s %12.6f" % ("Ratio edge/avg width",d[1]/avg_width))
    print("%20s %12.6f" % ("Ratio top/bottom",d[0]/d[2]))
    print("%20s %12.6f" % ("Ratio left/right",d[1]/d[3]))

    return d


def deg2m2(lon,lat):

    area_deg = (lat[1]-lat[0])*(lon[1]-lon[0])
    d = compute_distances(lon,lat)

    avg_width = (d[0] + d[2])/2.0
    avg_height = (d[1] + d[3])/2.0

    area_m2 = avg_width*avg_height

    return area_m2/area_deg


def region_coords(xll,xur, num_cells,lower,upper):
    """ Get location of region in coarse grid coordinates
    xll, xur     : lower and upper coordinates of approximate region
    num_clels    : Number of grid cells at coarsest level
    lower, upper :
    """

    # Get degrees per finest level :
    dx_coarse = (upper[0]-lower[0])/num_cells[0]
    dy_coarse = (upper[1]-lower[1])/num_cells[1]

    # Zoom region
    # xll = [-111.623926, 43.913661]
    # xur = [-111.620150, 43.916382]
    mx_xlow = np.floor((xll[0] - lower[0])/dx_coarse).astype(int)
    my_ylow = np.floor((xll[1] - lower[1])/dy_coarse).astype(int)
    xlow = lower[0] + mx_xlow*dx_coarse
    ylow = lower[1] + my_ylow*dy_coarse

    mx_zoom = np.ceil((xur[0] - xlow)/dx_coarse).astype(int)
    my_zoom = np.ceil((xur[1] - ylow)/dy_coarse).astype(int)

    region_lower = [xlow,ylow]
    region_upper = [xlow + mx_zoom*dx_coarse, ylow + my_zoom*dy_coarse]

    figsize = np.array([mx_zoom, my_zoom])   # [1,1]

    return region_lower, region_upper, figsize
