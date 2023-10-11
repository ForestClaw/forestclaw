/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
/** 
 * @file
 * @brief Functions needed for underlying wrapped types
 *        Most users will not need to include the file, unless doing something more advanced.
 */


#ifndef FCLAW2D_WRAP_H
#define FCLAW2D_WRAP_H

#include <fclaw_wrap.h>
#include <forestclaw2d.h>  /* Contains definition of patch-iterator callback */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

struct fclaw_domain* fclaw_domain_wrap_2d(fclaw2d_domain_t *domain2d);

void fclaw2d_patch_callback_wrap(fclaw2d_domain_t * domain, 
                                 fclaw2d_patch_t * patch,
                                 int blockno, int patchno, void *user);

void
fclaw2d_transfer_callback_wrap(fclaw2d_domain_t * old_domain,
                               fclaw2d_patch_t * old_patch,
                               fclaw2d_domain_t * new_domain,
                               fclaw2d_patch_t * new_patch,
                               int blockno,
                               int old_patchno, int new_patchno,
                               void *user);

void
fclaw2d_match_callback_wrap(fclaw2d_domain_t * old_domain,
                            fclaw2d_patch_t * old_patch_2d,
                            fclaw2d_domain_t * new_domain,
                            fclaw2d_patch_t * new_patch_2d,
                            fclaw2d_patch_relation_t newsize,
                            int blockno,
                            int old_patchno, int new_patchno,
                            void *user);

int 
fclaw2d_intersect_wrap (fclaw2d_domain_t * domain,
                        fclaw2d_patch_t * patch,
                        int blockno, int patchno,
                        void *ray, double *integral,
                        void *user);

int
fclaw2d_interpolate_point_wrap (fclaw2d_domain_t * domain,
                                fclaw2d_patch_t * patch,
                                int blockno, int patchno,
                                void *point, void *user);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
