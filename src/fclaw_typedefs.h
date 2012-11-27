#ifndef FCLAW_TYPEDEFS_H
#define FCLAW_TYPEDEFS_H

#include "fclaw_defs.H"
#include "amr_options.h"
#include "forestclaw2d.h"
#include "amr_mol.H"
#include "amr_single_step.H"

class ClawPatch;

typedef fclaw2d_level_time_data fclaw2d_level_time_data_t;

typedef void (*fclaw_level_advance_t)(fclaw2d_domain_t *domain,
                                         int level,
                                         fclaw2d_level_time_data_t *time_data);

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;
    double curr_time;

    // Always need a 'level advance' function - either
    fclaw_level_advance_t f_level_advance;

    // And either a single step solver or an mol solver.
    fclaw_single_step_patch_t f_single_step_patch;

    // MOL solver requires both a patch_rhs function and
    // an ODE solver (a mol_solver).
    fclaw_mol_rhs_patch_t f_mol_rhs_patch;
    fclaw_mol_solver_t f_mol_solver;
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
};


#endif
