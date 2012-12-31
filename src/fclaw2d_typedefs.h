#ifndef FCLAW2D_TYPEDEFS_H
#define FCLAW2D_TYPEDEFS_H

#include "fclaw2d_defs.H"
#include "amr_options.h"
#include "forestclaw2d.h"
#include "fclaw2d_solvers.H"

class ClawPatch;

typedef struct fclaw2d_level_time_data fclaw2d_level_time_data_t;

typedef void (*fclaw2d_problem_setup_t)(fclaw2d_domain_t* domain);

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;

    /* Some solver parms */
    void* waveprop_parms;

    /* Time at start of each subcycled time step */
    double curr_time;

    fclaw2d_problem_setup_t f_problem_setup;

    fclaw2d_solver_functions_t* solver_functions;

} fclaw2d_domain_data_t;

typedef struct fclaw2d_block_data
{
    int mthbc[FCLAW_NUMFACES];  /* >=0 for physical bc types */
} fclaw2d_block_data_t;

typedef struct fclaw2d_patch_data
{
    ClawPatch	*cp;
} fclaw2d_patch_data_t;

/* -----------------------------------------------------------
   Data needed for time stepping
   ----------------------------------------------------------- */
struct fclaw2d_level_time_data
{
    // Single step data. This always has to be set.
    double dt;
    double t_initial;
    double t_level;
    double t_coarse;

    // Needed for explicit CFL limited schemes
    double maxcfl;

    // Extra data that might be needed for more complicated time stepping
    // Not always set.
    double alpha;         // Fraction of coarser dt completed.
    double dt_coarse;
    bool is_coarsest;
    bool fixed_dt;
};


#endif
