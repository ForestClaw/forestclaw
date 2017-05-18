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

#ifndef FCLAW2D_DOMAIN_H
#define FCLAW2D_DOMAIN_H

#include <fclaw2d_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct fclaw2d_domain_data
{
    /* Debug counters and timers */
    int count_set_patch;
    int count_delete_patch;

    fclaw2d_domain_exchange_t *domain_exchange;
    fclaw2d_domain_indirect_t *domain_indirect;

} fclaw2d_domain_data_t;

void fclaw2d_domain_data_new(fclaw2d_domain_t *domain);

void fclaw2d_domain_data_delete(fclaw2d_domain_t* domain);

void fclaw2d_domain_setup(fclaw2d_global_t* glob,
                          fclaw2d_domain_t* new_domain);

void fclaw2d_domain_reset(fclaw2d_global_t* glob);

fclaw2d_domain_data_t*
fclaw2d_domain_get_data(fclaw2d_domain_t *domain);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif



#endif
