# forestclaw.py

from configparser import ConfigParser

class ForestClawData(object):
    # Several forestclaw attributes (ignore for now)

    def __init__(self, attributes=None):
        pass

    def write(self,rundata):
        geoclaw = ConfigParser(allow_no_value=True)

        clawdata = rundata.clawdata

        geoclaw['user'] = { 'example' : 0}
        geoclaw['clawpatch'] = {
            'mx' : clawdata.num_cells[0],
            'my' : clawdata.num_cells[1]}
        
        with open('geoclaw.ini','w') as geoclawfile:
            geoclaw.write(geoclawfile)
