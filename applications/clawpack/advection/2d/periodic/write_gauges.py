import fclaw_analysis

dim = 2;
gaugedata = fclaw_analysis.GaugeData(dim,min_time_increment=0)

# Periodic domain : [-1,1]x[-1,1]
#
# Format : [id, x, y, t0, t1]
#   id   : integer, identifying the gauge
#  x,y   : Location of the gauge
#  t0,t1 : (t0,t1) interval over which to monitor the gauge.

# This gauge travels in a circle and never leaves the domain
R = 0.5    # R < 1 : circle doesn't leave the domain. 
gaugedata.gauges.append([  0, R, 0,  0, 1.e10])

# This gauge travels in a straight line and may leave the domain eventually
# depending on the speed of the gauge
gaugedata.gauges.append([  1,  0,    0,    0, 1.e10])

# This gauge does not start in the domain, but may enter the domain 
# during the simulation. 
gaugedata.gauges.append([  2, -1.5, -1.5,  0, 1.e10])

# This is a static gauge - the "move" function doesn't move the gauge.  
gaugedata.gauges.append([20, -0.5, -0.5,  0, 1.e10])


gaugedata.write(data_source='write_gauges.py')
