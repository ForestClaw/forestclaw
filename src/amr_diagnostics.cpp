/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "amr_forestclaw.H"
#include "amr_utils.H"
#include "fclaw2d_convenience.h"
#include "fclaw_defs.H"

class ClawPatch;

// -----------------------------------------------------------------
// Diagnostics
// -----------------------------------------------------------------
static
void cb_check_conservation(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           void *user)
{
    ClawPatch *this_cp = get_clawpatch(this_patch);
    Real *sum = (Real*) user;

    *sum += this_cp->compute_sum();
}

void check_conservation(fclaw2d_domain_t *domain)
{
    Real sum = 0;
    fclaw2d_domain_iterate_patches(domain,cb_check_conservation,(void *) &sum);

    printf("Total sum = %24.16f\n",sum);
}

// Dump current patch
void cb_dump_patch(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patch = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patch + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump();
    }
}


void dump_patch(fclaw2d_domain_t *domain, int dump_patch)
{
    printf("Dumping patch (current) %d\n",dump_patch);
    fclaw2d_domain_iterate_patches(domain, cb_dump_patch, &dump_patch);
}

// Dump last patch
void cb_dump_last_patch(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patch = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patch + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump_last();
    }
}

void dump_last_patch(fclaw2d_domain_t *domain, int dump_patch)
{
    printf("Dumping patch (last) %d\n",dump_patch);
    fclaw2d_domain_iterate_patches(domain, cb_dump_last_patch,
                                   &dump_patch);
}

static
void cb_dump_time_interp_patch(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patch = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patch + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump_time_interp();
    }
}

void dump_time_interp_patch(fclaw2d_domain_t *domain, int dump_patch)
{
    printf("Dumping patch (time_interp) %d\n",dump_patch);
    fclaw2d_domain_iterate_patches(domain,
                                   cb_dump_time_interp_patch,
                                   &dump_patch);
}

static
void cb_dump_auxarray(fclaw2d_domain_t *domain,
	fclaw2d_patch_t *patch, int block_no, int patch_no, void *user)
{
    int dump_patchno = *((int *) user);
    int numb4 = domain->blocks[block_no].num_patches_before;
    if (patch_no == dump_patchno + numb4)
    {
        ClawPatch *cp = get_clawpatch(patch);
        cp->dump_auxarray();
    }
}

void dump_auxarray(fclaw2d_domain_t *domain, int dump_patchno)
{
    printf("Dumping patch (time_interp) %d\n",dump_patchno);
    fclaw2d_domain_iterate_patches(domain,
                                   cb_dump_auxarray,
                                   &dump_patchno);
}
