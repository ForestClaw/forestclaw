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

#ifndef FCLAW2D_PARTITION_H
#define FCLAW2D_PARTITION_H
#include <fclaw2d_vtable.h>

#if 0
#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif
#endif

void fclaw2d_partition_setup(fclaw2d_domain_t* domain);
void fclaw2d_partition_exchange_all(fclaw2d_domain_t* domain);

fclaw2d_domain_exchange_t *
fclaw2d_partition_get_exchange_data(fclaw2d_domain_t* domain);

/** Repartition all patches in parallel.
 * \param [in] mode             A level for amrinit, -1 for running simulation.
 */
void fclaw2d_partition_domain(fclaw2d_domain_t** domain, int mode);

void exchange_ghost_patch_data_levels(fclaw2d_domain_t* domain,
                                      int exchange_minlevel, int exchange_maxlevel);

void set_boundary_patch_ptrs(fclaw2d_domain_t* domain,int exchange_minlevel,
                             int exchange_maxlevel);

#if 0
#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
#endif
