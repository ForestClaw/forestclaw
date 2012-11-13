#ifndef FCLAW_TYPEDEFS_H
#define FCLAW_TYPEDEFS_H

#include "fclaw_defs.H"
#include "amr_options.h"
#include "forestclaw2d.h"

class ClawPatch;

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;
    double curr_time;
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


/*
typedef struct fclaw2d_time_interp_data
{
    int level;
    double t_level;
    double alpha;         // Fraction of coarser dt completed.
    double dt_fine;       // May require stages to get complete time step
    double dt_coarse;
    bool do_time_interp;
    bool is_coarsest;
} fclaw2d_time_interp_data_t;
*/

typedef struct fclaw2d_level_time_data
{
    double dt;
    double t;
    double dt_coarse;
    double maxcfl;
} fclaw2d_level_time_data_t;

typedef struct fclaw2d_patch_mol_data
{
    int count;
    int patch_size;
    fclaw2d_patch_t **patches; // ptrs to patches at a given level
    double *patch_data;      // Vectorized data
    double dx;
    double dy;
} fclaw2d_patch_mol_data_t;

typedef struct fclaw2d_f_exp_data_fort
{
    fclaw2d_domain_t* domain;
    int level;
    double t;
    double dt_coarse;
    fclaw2d_patch_mol_data_t *mol_data;
} fclaw2d_f_exp_data_fort_t;


#endif
