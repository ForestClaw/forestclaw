/*
  Copyright (c) 2019 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef FC2D_MULTIGRID_H
#define FC2D_MULTIGRID_H

#include <fclaw_base.h>   /* Needed for FCLAW_F77_FUNC */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw2d_patch;

typedef  struct fc2d_multigrid_vtable  fc2d_multigrid_vtable_t;


/* --------------------------- Fortran defs solver functions -------------------------- */

typedef  void (*fc2d_multigrid_fort_rhs_t)(const int* mbc,
                                           const int* mx, const int* my,
                                           const double* xlower, const double* ylower,
                                           const double* dx, const double* dy,
                                           double rhs[]);


/* -------------------------- Solver and utilities ------------------------------------ */

void fc2d_multigrid_solve(struct fclaw2d_global *glob);

/* --------------------------------- Virtual table ------------------------------------ */

struct fc2d_multigrid_vtable
{

	/* Fortran routines */
	fc2d_multigrid_fort_rhs_t       fort_rhs;	
	int is_set;

};

void fc2d_multigrid_solver_initialize(void);

fc2d_multigrid_vtable_t* fc2d_multigrid_vt(void);


/* ----------------------------- User access to solver functions ---------------------- */

void fc2d_multigrid_setprob(struct fclaw2d_global* glob);


void fc2d_multigrid_rhs(struct fclaw2d_global* glob,
                        struct fclaw2d_patch *patch,
                        int blockno,
                        int patchno);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_MULTIGRID_H */
