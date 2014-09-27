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
#include "fclaw2d_diagnostics.H"

static
void cb_check_conservation(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx,
                           void *user)
{
    ClawPatch *this_cp = get_clawpatch(this_patch);
    double *sum = (double*) user;

    *sum += this_cp->compute_sum();
}


void fclaw2d_check_sum(fclaw2d_domain_t *domain,
                       fclaw2d_cb_check_sum_t cb_check_sum)
{
    double sum = 0;
    fclaw2d_domain_iterate_patches(domain,cb_check_sum,(void *) &sum);

    sum = fclaw2d_domain_global_sum (domain, sum);
}


static
void cb_check_minimum(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    ClawPatch *this_cp = get_clawpatch(this_patch);
    double *minvalue = (double*) user;

    double mv = this_cp->compute_minimum();
    if (mv < *minvalue)
    {
        *minvalue = mv;
    }
}

void fclaw2d_check_minimum(fclaw2d_domain_t *domain,
                           fclaw2d_cb_check_minimum_t cb_check_minimum,
                           double *minvalue)
{
    /* minvalue should be set in calling function */
    fclaw2d_domain_iterate_patches(domain,cb_check_minimum,(void *) minvalue);

    *minvalue = fclaw2d_domain_global_minimum (domain, *minvalue);
}


static
void cb_check_maximum(fclaw2d_domain_t *domain,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx,
                      void *user)
{
    ClawPatch *this_cp = get_clawpatch(this_patch);
    double *maxvalue = (double*) user;

    double mv = this_cp->compute_max();
    if (mv > *maxvalue)
    {
        *maxvalue = mv;
    }
}


void fclaw2d_check_maximum(fclaw2d_domain_t *domain,
                           fclaw2d_cb_check_minimum_t cb_check_maximum,
                           double* maxvalue)
{
    /* Maxvalue should be set by the calling function */
    fclaw2d_domain_iterate_patches(domain,cb_check_maximum,(void *) maxvalue);

    *maxvalue = fclaw2d_domain_global_maximum (domain, *maxvalue);
}

void patch_diagnostics_dummy(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx, int this_patch_idx,
                             double t)
{
    /* Don't do anything */
}



void initialize_diagnostic_functions(fclaw2d_diagnostic_functions_t* diagnostic_functions)
{
    diagnostic_functions->f_patch_diagnostics = patch_diagnostics_dummy;
}
