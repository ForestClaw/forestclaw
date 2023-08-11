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

#ifndef FCLAW2D_DOMAIN_H
#define FCLAW2D_DOMAIN_H

#include <forestclaw2d.h>  /* Needed for domain_exchange/domain_indirect info */
#include <fclaw_domain.h>

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
struct fclaw2d_global;

typedef struct fclaw2d_domain_data
{
    /* Debug counters and timers */
    int count_set_patch;
    int count_delete_patch;

    fclaw2d_domain_exchange_t *domain_exchange;
    fclaw2d_domain_indirect_t *domain_indirect;

} fclaw2d_domain_data_t;

void fclaw2d_domain_data_new(struct fclaw2d_domain *domain);

void fclaw2d_domain_data_delete(struct fclaw2d_domain* domain);

void fclaw2d_domain_setup(struct fclaw2d_global* glob,
                          struct fclaw2d_domain* new_domain);

void fclaw2d_domain_reset(struct fclaw2d_global* glob);

fclaw2d_domain_data_t* fclaw2d_domain_get_data(struct fclaw2d_domain *domain);

/* OpenMP iterator (not part of forestclaw2d.h */
void fclaw2d_domain_iterate_level_mthread (struct fclaw2d_domain * domain, int level,
                                           fclaw2d_patch_callback_t pcb, void *user);

/* below are the functions needed for dimension independence */

/** safeguard value for dimension-independent domain */
#define FCLAW2D_DOMAIN_MAGIC 0x56780202

void fclaw2d_domain_iterate_cb
  (fclaw2d_domain_t * d2, fclaw_patch_t * patch,
   int blockno, int patchno, void *user);

/** Construct a dimension-specific domain and initialize patch data.
 * \param [in] domain           Freshly created, valid domain structure
 *                              without any user data attached to it.
 * \param [in] init             This callback is pointed to a freshly
 *                              allocated fclaw_patch whose 2D data it
 *                              is supposed to fill with fitting values.
 *                              This includes allocating and filling
 *                              the patch user pointer, which the
 *                              callback can do using the third
 *                              parameter to this function, \a user.
 *                              May be NULL for no callback.
 * \param [in,out] user         Pointer passed through to \a init.
 * \return                      Initialized dimension-independent domain.
 */
fclaw_domain_t *fclaw_domain_new2d (fclaw2d_domain_t * domain,
                                    fclaw_domain_callback_t init, void *user);

/** Destruct a dimension-specific domain and deinitialize patch data.
 * \param [in] domain           Valid domain-independent domain structure.
 * \param [in] dele             This callback is pointed to each existing
 *                              fclaw_patch whose 2D data it is supposed to
 *                              deinitialize, as well as the patch user data.
 *                              Note, just don't allocate the patch itself!
 *                              May be NULL for no callback.
 * \param [in,out] user         Pointer passed through to \a init.
 */
void fclaw_domain_destroy2d (fclaw_domain_t * domain,
                             fclaw_domain_callback_t dele, void *user);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* FCLAW2D_DOMAIN_H */
