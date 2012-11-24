#ifndef FCLAW_TYPEDEFS_H
#define FCLAW_TYPEDEFS_H

#include "fclaw_defs.H"
#include "amr_options.h"
#include "forestclaw2d.h"
#include "amr_mol.H"

class ClawPatch;

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;
    double curr_time;
    fclaw_mol_solver_t f_mol_solver;
    fclaw_mol_rhs_patch_t f_mol_rhs_patch;
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
typedef
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
} fclaw2d_level_time_data_t;


/*
typedef struct fclaw2d_methods
{
    // static void* f_single_step_update = NULL;
    fclaw_mol_solver_t *f_mol_solver;
    fclaw_mol_rhs_patch_t *f_mol_rhs_patch;
} fclaw2d_methods_t;
*/

#endif
