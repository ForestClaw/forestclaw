import fclaw_analysis

dim = 2;
gaugedata = fclaw_analysis.GaugeData(dim,min_time_increment=0)

# Periodic domain : [-1,1]x[-1,1]
#
# Format : [id, x, y, t0, t1]
#   id   : integer, identifying the gauge
#  x,y   : Location of the gauge
#  t0,t1 : (t0,t1) interval over which to monitor the gauge.

gaugedata.gauges.append([  0,  0,    0,    0, 1.e10])
gaugedata.gauges.append([  1, -0.5, -0.5,  0, 1.e10])
gaugedata.gauges.append([  2,  0.5,  0.5,  0, 1.e10])

# A gauge that is not in the domain (for testing purposes)
gaugedata.gauges.append([  4, 0, -2,  0, 1.e10])

gaugedata.write(data_source='write_gauges.py')
