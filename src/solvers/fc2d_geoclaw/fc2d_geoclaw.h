/*
  Copyright (c) 2012 Carsten Burstedde, Donna Calhoun, Yu-Hsuan Shih
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

#ifndef FC2D_GEOCLAW_H
#define FC2D_GEOCLAW_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */
#include <fc2d_geoclaw_fort.h>  /* Needed for virtual functions */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fc2d_geoclaw_vtable fc2d_geoclaw_vtable_t;

/* Forward declarations */
struct fclaw2d_patch_transform_data;
struct fclaw2d_global;
struct fclaw2d_patch;
struct geoclaw_gauge;


void fc2d_geoclaw_run(struct fclaw2d_global *glob);

/* ------------------------------------- Access functions ---------------------------------- */

void fc2d_geoclaw_module_setup(struct fclaw2d_global *glob);

void fc2d_geoclaw_output(struct fclaw2d_global *glob, int iframe);

/* ------------------------------------- Virtual table ----------------------------------- */

void fc2d_geoclaw_solver_initialize();

fc2d_geoclaw_vtable_t* fc2d_geoclaw_vt();


struct fc2d_geoclaw_vtable
{
    fc2d_geoclaw_setprob_t  setprob;
    fc2d_geoclaw_bc2_t      bc2;
    fc2d_geoclaw_qinit_t    qinit;
    fc2d_geoclaw_setaux_t   setaux;
    fc2d_geoclaw_b4step2_t  b4step2;
    fc2d_geoclaw_src2_t     src2;
    fc2d_geoclaw_rpn2_t     rpn2;
    fc2d_geoclaw_rpt2_t     rpt2;
    fc2d_geoclaw_fluxfun_t  fluxfun;

    int is_set;
};


#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_CLAWPACH5_H */
