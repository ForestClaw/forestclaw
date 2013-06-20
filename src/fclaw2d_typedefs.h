#ifndef FCLAW2D_TYPEDEFS_H
#define FCLAW2D_TYPEDEFS_H

#include "fclaw2d_defs.H"
#include "amr_options.h"
#include "forestclaw2d.h"
#include "fclaw2d_solvers.H"

class ClawPatch;

typedef struct fclaw2d_level_time_data fclaw2d_level_time_data_t;

typedef void (*fclaw2d_problem_setup_t)(fclaw2d_domain_t* domain);

typedef void (*fclaw2d_patch_output_t)(fclaw2d_domain_t* domain, fclaw2d_patch_t *this_patch,
                                       int this_block_idx, int this_patch_idx,
                                       int iframe, int num, int matlab_level);

typedef fclaw_bool (*fclaw2d_patch_tag4refinement_t)(fclaw2d_domain_t *domain,
                                                     fclaw2d_patch_t *this_patch,
                                                     int this_block_idx, int this_patch_idx,
                                                     int initflag);

typedef fclaw_bool (*fclaw2d_patch_tag4coarsening_t)(fclaw2d_domain_t *domain,
                                                     fclaw2d_patch_t *sibling_patch,
                                                     int this_block_idx,
                                                     int sibling0_patch_idx,
                                                     ClawPatch *cp_new_coarse);


typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;

    /* Some solver parms */
    void* waveprop_parms;

    void* manyclaw_parms;

    /* Time at start of each subcycled time step */
    double curr_time;

    fclaw2d_problem_setup_t f_problem_setup;
    fclaw2d_patch_output_t f_patch_output;
    fclaw2d_patch_tag4refinement_t f_patch_tag4refinement;
    fclaw2d_patch_tag4coarsening_t f_patch_tag4coarsening;

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
