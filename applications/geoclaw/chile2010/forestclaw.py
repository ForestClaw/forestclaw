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
        geoclaw.optionxform = str # To original case of comments. 

        clawdata = rundata.clawdata
        geo_data = rundata.geo_data
        refinement_data = rundata.refinement_data
        amrdata = rundata.amrdata

        geoclaw['user'] = {
        '   example' : 0
        }
        
        geoclaw['clawpatch'] = {
        '# Grid dimensions' : None,
        '   # mx_leaf' : None,
        '   mx' : clawdata.num_cells[0],"\n"
        '   # my_leaf' :None, 
        '   my' : clawdata.num_cells[1], "\n"
        '   # Number of ghost cells': None,
        '   mbc': clawdata.num_ghost, "\n"

        '   maux': clawdata.num_aux,"\n"
        '   # Number of equations' : None,
        '   meqn': clawdata.num_eqn        
            }

        geoclaw['Options'] = {

        '# Regridding information' : None,
        '   # Minimum level':None,
        '   minlevel' : self.minlevel,"\n"
        '   # Maximum levels of refinement':None,
        '   maxlevel' : self.maxlevel,"\n"
        "  # Regrid every 'regrid_interval' time steps.":None,
        '   weighted_partition' : self.weighted_partition,"\n"

        '   regrid-interval' : self.regrid_interval,
        '   refine-threshold' : self.refine_threshold,
        '   coarsen-threshold' : self.coarsen_threshold,
        '   smooth-refine' : self.smooth_refine,
        '   smooth-level' : self.smooth_level,
        '   outstyle-uses-maxlevel' : self.outstyle_uses_maxlevel,
        '   subcycle' : self.subcycle,
        '   advance-one-step' : self.advance_one_step,"\n"

        '   # Diagnostics and debugging':None,
        '   # Trap floating point errors.':None,
        '   trapfpe' : self.trapfpe,"\n"
        '   # Attach mpi processes in gdb':None,
        '   mpi-debug' : self.mpi_debug,"\n"
        '   conservation-check' : self.conservation_check,
        '   run-user-diagnostics' : self.run_user_diagnostics,"\n"

        '# File and console IO' : None,
        '   output-gauges' : self.output_gauges,
        '   output' : self.output,
        '   verbosity' : self.verbosity,"\n"

        '# Time stepping' : None, 
        '   # Final time':None,
        '   tfinal': clawdata.tfinal, "\n"     

        '   # Take a fixed time step':None,
        '   use-fixed-dt': not clawdata.dt_variable,"\n"
        "  # Initial time step for 'minlevel'":None,
        '   initial-dt': clawdata.dt_initial,"\n"
        '   # maximum cfl':None,
        '   max-cfl': clawdata.cfl_max,"\n"
        '   # desired cfl':None,
        '   desired-cfl': clawdata.cfl_desired,"\n"

        '   # 1 : Output steps  = tfinal/nout;':None,                  
        '   # 2 : not implemented;':None,                              
        '   # 3 : Take nout steps;  save files every nstep steps.':None,
        '   outstyle': clawdata.output_style,"\n" 

        '   # Used for all three out styles;  has different meaning, though':None,                                     
        '   nout': clawdata.num_output_times,"\n"
        '   # Only used if outstyle is 3':None,
        '   nstep': clawdata.output_style,"\n"

        '# Mapping' : None,"\n"
        '   # Domain dimensions' : None,
        '   ax': clawdata.lower[0],
        '   bx': clawdata.upper[0],
        '   ay': clawdata.lower[1],
        '   by': clawdata.upper[1]
        }
        
        #mthbc
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

        #Apply same idea for the limiter! 
        lim_in  = [ s for s in clawdata.limiter]

        lim = [0]*clawdata.num_waves
        lim_str = ""
        bc_dict = {'none' : 0, 'minmod' : 1, 'superbee' : 2, 'vanleer' : 3, 'mc' : 4}
        for k in range(clawdata.num_waves):
            if type(lim_in[k]) == str:
                lim[k] = bc_dict[lim_in[k]]
            else:
                lim[k] = lim_in[k]

            lim_str += " " + str(lim[k])

        # Apply same idea to order 
        ord_in  = clawdata.order

        ord = [0]*clawdata.order
        ord_str = ""
        bc_dict = {'godunov' : 1, 'Lax-Wendroff' : 2}
        for k in range(clawdata.order):
            if type(ord_in) == str:
                ord[k] = bc_dict[ord_in]
            else:
                ord[k] = ord_in

            ord_str += " " + str(ord[k])

        # Can we get dashes instead of underscores?   yes


        geoclaw['geoclaw'] = {
            '   # normal and transverse order': None,
            '   # Order of accuracy:': None,
            '     # 1 => Godunov,': None,  
            '     # 2 => Lax-Wendroff plus limiters': None,
            #order': "2 2", "\n" 
            '   order': ord_str,"\n"

            '   # mcapa' : None,
            '   mcapa': clawdata.capa_index,"\n"
            '   # mbathy' : None,
            '   mbathy': 1,"\n"
            '   # src_term' : None,
            '   src-term': 1,"\n"

            '   # mwaves': None,
            '   mwaves': clawdata.num_waves,"\n"

            "  # mthlim (is a vector in general, with 'mwaves' entries": None,
            '    # List of limiters to use for each wave family:': None,
            '    # Required:  len(limiter) == num_waves': None,
            '    # Some options:': None,
            "    #   0 or 'none'     ==> no limiter (Lax-Wendroff)": None,
            "    #   1 or 'minmod'   ==> minmod": None,
            "    #   2 or 'superbee' ==> superbee": None,
            "    #   3 or 'mc'       ==> MC limiter": None,
            "    #   4 or 'vanleer'  ==> van Leer": None,
            #'   mthlim': "4 4 4","\n"
            '   mthlim' : lim_str,"\n"

            '   # mthbc (=left,right,bottom,top)' : None,
                # 'mthbc': "1  1  1  1",   # 0,1,2,3
            '   # Choice of BCs at xlower and xupper:':None,
            '   # 0 => user specified (must modify bcN.f to use this option)':None,
            '   # 1 => extrapolation (non-reflecting outflow)':None,
            '   # 2 => periodic (must specify this at both boundaries)':None,
            '   # 3 => solid wall for systems where q(2) is normal velocity':None,
            '   mthbc' : mthbc_str,"\n"

            '# Coarsening' : None,            #(probably going to go away)
            '   dry-tolerance-c': geo_data.dry_tolerance,
            '   wave-tolerance-c': refinement_data.wave_tolerance,
            '   speed-tolerance-entries-c': 6,
            '   speed-tolerance-c': "1e12 1e12 1e12 1e12 1e12 1e12","\n"

            '   # Output' : None,
            '   ascii-out': True,"\n"

            '   # Number of lines in gauge file to store in memory before printing' : None,
            '   gauge-buffer-length': 100

            }
            
        with open('geoclaw.ini','w') as geoclawfile:
            geoclaw.write(geoclawfile)
