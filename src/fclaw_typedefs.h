#ifndef FCLAW_TYPEDEFS_H
#define FCLAW_TYPEDEFS_H

#include "fclaw_defs.H"

#include "forestclaw2d.h"
#include "fclaw2d_convenience.h"
#include "amr_options.h"

class ClawPatch;

typedef struct fclaw2d_domain_data
{
    const amr_options_t *amropts;
    double curr_time;
} fclaw2d_domain_data_t;

typedef struct fclaw2d_block_data
{
    int mthbc[4];               /* >=0 for physical bc types */
} fclaw2d_block_data_t;

typedef struct fclaw2d_patch_data
{
    ClawPatch	*cp;
} fclaw2d_patch_data_t;

typedef struct fclaw2d_level_time_data
{
    double dt;
    double maxcfl;
    double t;
} fclaw2d_level_time_data_t;

#endif
