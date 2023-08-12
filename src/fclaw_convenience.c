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

#include "fclaw2d_convenience.h"
#include "fclaw3d_convenience.h"

void fclaw2d_domain_destroy (fclaw_domain_t * domain);

fclaw_domain_t *fclaw2d_domain_adapt (fclaw_domain_t * domain);

fclaw_domain_t *fclaw2d_domain_partition (fclaw_domain_t * domain,
                                            int weight_exponent);

void fclaw2d_domain_partition_unchanged (fclaw_domain_t * domain,
                                         int *unchanged_first,
                                         int *unchanged_length,
                                         int *unchanged_old_first);

void fclaw2d_domain_complete (fclaw_domain_t * domain);

void fclaw2d_domain_write_vtk (fclaw_domain_t * domain,
                               const char *basename);

void fclaw2d_domain_list_levels (fclaw_domain_t * domain, int log_priority);

void fclaw2d_domain_list_neighbors (fclaw_domain_t * domain,
                                    int log_priority);

void fclaw2d_domain_list_adapted (fclaw_domain_t * old_domain,
                                  fclaw_domain_t * new_domain,
                                  int log_priority);


void fclaw3d_domain_destroy (fclaw_domain_t * domain);

fclaw_domain_t *fclaw3d_domain_adapt (fclaw_domain_t * domain);

fclaw_domain_t *fclaw3d_domain_partition (fclaw_domain_t * domain,
                                            int weight_exponent);

void fclaw3d_domain_partition_unchanged (fclaw_domain_t * domain,
                                         int *unchanged_first,
                                         int *unchanged_length,
                                         int *unchanged_old_first);

void fclaw3d_domain_complete (fclaw_domain_t * domain);

void fclaw3d_domain_write_vtk (fclaw_domain_t * domain,
                               const char *basename);

void fclaw3d_domain_list_levels (fclaw_domain_t * domain, int log_priority);

void fclaw3d_domain_list_neighbors (fclaw_domain_t * domain,
                                    int log_priority);

void fclaw3d_domain_list_adapted (fclaw_domain_t * old_domain,
                                  fclaw_domain_t * new_domain,
                                  int log_priority);


void fclaw_domain_destroy (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_destroy(domain);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_destroy(domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

fclaw_domain_t* fclaw_domain_adapt (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        return fclaw2d_domain_adapt(domain);
    }
    else if(domain->dim == 3)
    {
        return fclaw3d_domain_adapt(domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

fclaw_domain_t* fclaw_domain_partition (fclaw_domain_t * domain,
                                            int weight_exponent)
{
    if(domain->dim == 2)
    {
        return fclaw2d_domain_partition(domain,weight_exponent);
    }
    else if(domain->dim == 3)
    {
        return fclaw3d_domain_partition(domain,weight_exponent);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_partition_unchanged (fclaw_domain_t * domain,
                                         int *unchanged_first,
                                         int *unchanged_length,
                                         int *unchanged_old_first)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_partition_unchanged(domain,unchanged_first,
                                           unchanged_length,
                                           unchanged_old_first);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_partition_unchanged(domain,unchanged_first,
                                           unchanged_length,
                                           unchanged_old_first);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_complete (fclaw_domain_t * domain)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_complete(domain);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_complete(domain);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_write_vtk (fclaw_domain_t * domain,
                             const char *basename)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_write_vtk(domain,basename);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_write_vtk(domain,basename);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_levels (fclaw_domain_t * domain, int log_priority)
{
    if(domain->dim == 2)
    {
        fclaw_domain_list_levels(domain,log_priority);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_list_levels(domain,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_neighbors (fclaw_domain_t * domain, int log_priority)
{
    if(domain->dim == 2)
    {
        fclaw_domain_list_neighbors(domain,log_priority);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_list_neighbors(domain,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_adapted (fclaw_domain_t * old_domain,
                                fclaw_domain_t * new_domain,
                                int log_priority)
{
    if(old_domain->dim == 2)
    {
        fclaw2d_domain_list_adapted(old_domain,new_domain,log_priority);
    }
    else if(old_domain->dim == 3)
    {
        fclaw3d_domain_list_adapted(old_domain,new_domain,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}