# forestclaw.py

from configparser import ConfigParser

class ForestClawData(object):
    # Several forestclaw attributes (ignore for now)

    def __init__(self):
        self.minlevel = 0
        self.maxlevel = 0

        self.weighted_partition = True
        self.regrid_interval = 1
        self.refine_threshold = 0.5
        self.coarsen_threshold = 0.1
        self.smooth_refine = False
        self.smooth_level = 0
        self.outstyle_uses_maxlevel = False
        self.subcycle = False
        self.advance_one_step = False
        self.trapfpe = False
        self.mpi_debug = False
        self.conservation_check = False
        self.run_user_diagnostics = False
        self.output_gauges = False
        self.output = True
        self.verbosity = 'essential'
        self.report_timing=True
        self.report_timing_verbosity = 'wall'


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

            'minlevel' : self.minlevel,
            'maxlevel' : self.maxlevel,
            'weighted_partition' : self.weighted_partition,
            'regrid_interval' : self.regrid_interval,
            'refine_threshold' : self.refine_threshold,
            'coarsen_threshold' : self.coarsen_threshold,
            'smooth_refine' : self.smooth_refine,
            'smooth_level' : self.smooth_level,
            'outstyle_uses_maxlevel' : self.outstyle_uses_maxlevel,
            'subcycle' : self.subcycle,
            'advance_one_step' : self.advance_one_step,
            'trapfpe' : self.trapfpe,
            'mpi_debug' : self.mpi_debug,
            'conservation_check' : self.conservation_check,
            'run_user_diagnostics' : self.run_user_diagnostics,
            'output_gauges' : self.output_gauges,
            'output' : self.output,
            'verbosity' : self.verbosity,
            #Time stepping 
            'tfinal': clawdata.tfinal,      # Final time
            'use_fixed_dt': not clawdata.dt_variable,
            'initial_dt': clawdata.dt_initial,
            'max_cfl': clawdata.cfl_max,
            'desired_cfl': clawdata.cfl_desired,
            'outstyle': clawdata.output_style,
            'nout': clawdata.num_output_times,
            'nstep': clawdata.output_step_interval,
            # Mapping
            # Domain dimensions
            'ax': clawdata.lower[0],
            'bx': clawdata.upper[0],
            'ay': clawdata.lower[1],
            'by': clawdata.upper[1]
            }


# Choice of BCs at xlower and xupper:
#   0 => user specified (must modify bcN.f to use this option)
#   1 => extrapolation (non-reflecting outflow)
#   2 => periodic (must specify this at both boundaries)
#   3 => solid wall for systems where q(2) is normal velocity

#clawdata.bc_lower[0] = 'extrap'
#clawdata.bc_upper[0] = 'extrap'
#clawdata.bc_lower[1] = 'extrap'
#clawdata.bc_upper[1] = 'extrap'

        
        mthbc_in  = [clawdata.bc_lower[0], clawdata.bc_upper[0], 
                     clawdata.bc_lower[1], clawdata.bc_upper[1]]
        mthbc = [0]*4
        mthbc_str = ""
        bc_dict = {'user' : 0, 'extrap' : 1, 'periodic' : 2, 'wall' : 3}
        for k in range(4):
            if type(mthbc_in[k]) == str:
                mthbc[k] = bc_dict[mthbc_in[k]]
            else:
                mthbc[k] = mthbc_in[k]

            mthbc_str += " " + str(mthbc[k])

        # Apply same idea for the limiter! 

        # Apply same idea to order 

        # Can we get dashes instead of underscores?    


        geoclaw['geoclaw'] = {
                'order': "2 2   # Order of the method",  
                'mcapa': clawdata.capa_index,
                'mbathy': 1,
                'src_term': 1,
                'mwaves': clawdata.num_waves,
                'mthlim': "4 4 4",
                # 'mthbc': "1  1  1  1",   # 0,1,2,3
                'mthbc' : mthbc_str,

                # Coarsening (probably going to go away)
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
