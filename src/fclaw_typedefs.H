#ifndef FCLAW_TYPEDEFS_H
#define FCLAW_TYPEDEFS_H

#include "fclaw2d_capi.h"
#include "fclaw2d_single_step.h"

#include "fclaw_defs.H"
#include "amr_mol.H"

class ClawPatch;

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;
    double curr_time;

    // Always need a 'level advance' function - either
    fclaw2d_level_advance_t f_level_advance;

    // And either a single step solver or an mol solver.
    fclaw2d_single_step_patch_t f_single_step_patch;

    // MOL solver requires both a patch_rhs function and
    // an ODE solver (a mol_solver).
    fclaw_mol_rhs_patch_t f_mol_rhs_patch;
    fclaw_mol_solver_t f_mol_solver;
}
fclaw2d_domain_data_t;

typedef struct fclaw2d_block_data
{
    int mthbc[FCLAW_NUMFACES];  /* >=0 for physical bc types */
}
fclaw2d_block_data_t;

typedef struct fclaw2d_patch_data
{
    ClawPatch	*cp;
}
fclaw2d_patch_data_t;

#endif
