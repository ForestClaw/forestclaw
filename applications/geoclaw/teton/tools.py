import numpy as np

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
