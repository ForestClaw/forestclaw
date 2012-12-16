#ifndef FCLAW2D_TYPEDEFS_H
#define FCLAW2D_TYPEDEFS_H

#include "fclaw2d_defs.H"
#include "amr_options.h"
#include "forestclaw2d.h"
#include "amr_solver_typedefs.H"

class ClawPatch;

typedef fclaw2d_level_time_data fclaw2d_level_time_data_t;

typedef void (*fclaw2d_level_advance_t)(fclaw2d_domain_t *domain,
                                        int level,
                                        fclaw2d_level_time_data_t *time_data);

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;

    /* Time at start of each subcycled time step */
    double curr_time;

    fclaw2d_patch_setup_t f_patch_setup_ptr;
    fclaw2d_patch_initialize_t f_patch_initialize_ptr;

    /* A single step solver.  Note that the user may want to
       change this one out */
    fclaw2d_single_step_level_t f_single_step_level_ptr;
    fclaw2d_single_step_update_patch_t f_single_step_update_patch_ptr;

    /* MOL solver requires both a patch_rhs function and
       an ODE solver (a mol_solver).  Note that the
       function 'fclaw2d_mol_step_level' is fixed;  I don't
       expect users to change this (unlike with
       'fclaw2d_single_step_level_t */

    fclaw2d_ode_solver_level_t f_ode_solver_level_ptr;
    fclaw2d_ode_solver_rhs_patch_t f_ode_solver_rhs_patch_ptr;
    fclaw2d_patch_physbc_t f_patch_physbc_ptr;

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
