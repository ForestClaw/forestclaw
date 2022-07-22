/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef ADVECTION_USER_H
#define ADVECTION_USER_H

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch_pillow.h>

/* Headers for both Clawpack 4.6 and  Clawpack 5.0 */
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_clawpatch_fort.h>


/* Clawpack 4.6 headers */
#include <fc2d_clawpack46.h>  
#include <fc2d_clawpack46_options.h>
#include <fc2d_clawpack46_fort.h>  
#include <clawpack46_user_fort.h>  
#include <fclaw2d_clawpatch46_fort.h>


/* Clawpack 5.0 headers */
#include <fc2d_clawpack5.h>
#include <fc2d_clawpack5_options.h>
#include <fc2d_clawpack5_fort.h>
#include <clawpack5_user_fort.h>
#include <fclaw2d_clawpatch5_fort.h>


/* Headers for common FORTRAN files */
#include "advection_user_fort.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* Common patch setup - this calls a common SETAUX routine that sets 
   velocities based on a streamfunction */

void advection_patch_setup_manifold(fclaw2d_global_t *glob,
                                    fclaw2d_patch_t *patch,
                                    int block,
                                    int patchno,
                                    int claw_version);

void advection_b4step2_manifold(fclaw2d_global_t *glob,
                                fclaw2d_patch_t *patch,
                                int blockno,
                                int patchno,
                                double t,
                                double dt,
                                int claw_version);

/* -------------------------------------- No map -------------------------------------- */

fclaw2d_map_context_t* fclaw2d_map_new_nomap();


/* --------------------------------- Square mappings ---------------------------------- */

fclaw2d_map_context_t* fclaw2d_map_new_identity(fclaw2d_map_context_t *brick);

fclaw2d_map_context_t* fclaw2d_map_new_cart(fclaw2d_map_context_t* brick,
                                            const double scale[],
                                            const double shift[]);
  
fclaw2d_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double alpha);

fclaw2d_map_context_t* fclaw2d_map_new_squareddisk(const double scale[],
                                                   const double shift[],
                                                   const double alpha);
  
fclaw2d_map_context_t* fclaw2d_map_new_bilinear(fclaw2d_map_context_t *brick,
                                                const double scale[],
                                                const double shift[],
                                                const double center[]);

/* ---------------------------------- Disk mappings ----------------------------------- */

fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk(const double scale[],
                                                  const double shift[],
                                                  const double rotate[]);

fclaw2d_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);

/* --------------------------------- Annulus mapping ---------------------------------- */
fclaw2d_map_context_t *
    fclaw2d_map_new_annulus (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double rotate[],
                             const double alpha,
                             const double theta[]);


/* --------------------------------- Latlong mapping ---------------------------------- */
fclaw2d_map_context_t *
    fclaw2d_map_new_latlong (fclaw2d_map_context_t* brick,
                             const double scale[],
                             const double lat[],
                             const double longitude[],
                             const int a, const int b);

/* --------------------------------- Hemisphere mappings ------------------------------ */

fclaw2d_map_context_t* fclaw2d_map_new_pillowsphere5(const double scale[],
                                                     const double rotate[],
                                                     const double alpha);

/* --------------------------------- Sphere mappings ---------------------------------- */

fclaw2d_map_context_t* fclaw2d_map_new_pillowsphere(const double scale[],
                                                    const double rotate[]);

fclaw2d_map_context_t * fclaw2d_map_new_cubedsphere (const double scale[],
                                                     const double rotate[]);

fclaw2d_map_context_t * fclaw2d_map_new_pillowsphere (const double scale[],
                                                      const double rotate[]);


/* --------------------------------- Torus mappings ---------------------------------- */
fclaw2d_map_context_t *
    fclaw2d_map_new_torus (fclaw2d_map_context_t* brick,
                           const double scale[],
                           const double rotate[],
                           const double alpha,
                           const double beta);








#ifdef __cplusplus
}
#endif

#endif
