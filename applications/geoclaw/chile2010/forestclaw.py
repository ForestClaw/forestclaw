# forestclaw.py

from configparser import ConfigParser

class ForestClawData(object):
    # Several forestclaw attributes (ignore for now)

    def __init__(self, attributes=None):
        pass

    def write(self,rundata):
        geoclaw = ConfigParser(allow_no_value=True)

        clawdata = rundata.clawdata
        geo_data = rundata.geo_data
        refinement_data = rundata.refinement_data
        amrdata = rundata.amrdata

        geoclaw['user'] = { 'example' : 0}
        geoclaw['clawpatch'] = {
            #Grid dimensions
            'mx' : clawdata.num_cells[0],   #mx_leaf
            'my' : clawdata.num_cells[1],   #my_leaf
            'mbc': clawdata.num_ghost,      #Number of ghost cells
            'maux': clawdata.num_aux,
            'meqn': clawdata.num_eqn        #Number of equations
            }

        geoclaw['Options'] = {
            'minlevel' : 0,
            'maxlevel' : refinement_data.max_level_deep,  #Maximmum level of refinement
            'regrid_interval': amrdata.regrid_interval,
            'refine_threshold': 0.005,
            'coarsen_threshold': False,
            'smooth-level': 4,
            'coarsen-delay': 0,
            #Time stepping 
            'tfinal': clawdata.tfinal,      # Final time
            'use_fixed_dt': clawdata.dt_variable,
            'initial_dt': clawdata.dt_initial,
            'max_cfl': clawdata.cfl_max,
            'desired_cfl': clawdata.cfl_desired,
            'outstyle': clawdata.output_style,
            'nout': clawdata.num_output_times,
            'nstep': clawdata.output_step_interval,
            'subcycle': False,
            'outstyle-uses-maxlevel': True,
            'weighted_partition': True,
            'advanced-one-step': True,
            # File and console IO
            'verbosity': clawdata.verbosity,
            'output-gauges': True,
            'output': clawdata.output_t0 ,
            #Diagnostics and dubbing
            'trapfpe': False,
            'mpi_debug': False,
            'conservation-check': False,
            'run-user-diagnostics': False,
            # Mapping
            # Domain dimensions
            'ax': clawdata.lower[0],
            'bx': clawdata.upper[0],
            'ay': clawdata.lower[1],
            'by': clawdata.upper[1]
            }

        geoclaw['geoclaw'] = {
                'order': [2,2],
                'mcapa': clawdata.capa_index,
                'mbathy': 1,
                'src_term': 1,
                'mwaves': clawdata.num_waves,
                'mthlim': [4, 4, 4],
                'mthbc': [1, 1, 1, 1],
                # Coarsening
                'dry_tolerance_c': geo_data.dry_tolerance,
                'wave_tolerance_c': refinement_data.wave_tolerance,
                'speed_tolerance_entries_c': 6,
                'speed_tolerance_c': [1e12, 1e12, 1e12, 1e12, 1e12, 1e12],
                # Output
                'ascii-out': True,
                #Number of lines in gauge file to store in memory before printing
                'gauge-buffer-length': 100
            }
            
        with open('geoclaw.ini','w') as geoclawfile:
            geoclaw.write(geoclawfile)
