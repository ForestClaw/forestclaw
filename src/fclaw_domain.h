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
/** 
 * @file
 * Domain structures and routines
 */

#ifndef FCLAW_DOMAIN_H
#define FCLAW_DOMAIN_H

#include <forestclaw.h>  /* Needed for domain_exchange/domain_indirect info */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* CB: chicken and egg -- should global include domain or vice versa?
       Believe removing any dependence on global from domain will work.
       Setting a global timer in domain_setup may likely be refactored.
       Deleting patch and exchange data in domain_reset might go into a
       toplevel algorithmic function quite naturally outside of domain.
 */
struct fclaw_global;

void fclaw_domain_setup(struct fclaw_global* glob,
                          struct fclaw_domain* new_domain);

void fclaw_domain_reset(struct fclaw_global* glob);

/* OpenMP iterator (not part of forestclaw2d.h */
void fclaw_domain_iterate_level_mthread (struct fclaw_domain * domain, int level,
                                           fclaw_patch_callback_t pcb, void *user);

/**
 * @brief Check if domain has exchange data allocated
 * 
 * @param domain the domain
 * @return int true if exchange data allocated
 */
int fclaw_domain_exchange_allocated(fclaw_domain_t *domain);

/* below are the functions needed for dimension independence */

/** safeguard value for dimension-independent domain */
#define FCLAW2D_DOMAIN_MAGIC 0x56780202

void fclaw2d_domain_iterate_cb
  (fclaw_domain_t * d2, fclaw_patch_t * patch,
   int blockno, int patchno, void *user);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* FCLAW2D_DOMAIN_H */
