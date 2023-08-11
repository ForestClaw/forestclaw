/*
  Copyright (c) 2019-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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

#ifndef FC2D_THUNDEREGG_H
#define FC2D_THUNDEREGG_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;
struct fclaw_patch;

typedef  struct fc2d_thunderegg_vtable  fc2d_thunderegg_vtable_t;



/* --------------------------- Fortran defs solver functions -------------------------- */

typedef  void (*fc2d_thunderegg_fort_rhs_t)(const int* blockno, 
                                           const int* mbc,
                                           const int* mx, const int* my,
                                           const int* mfields,
                                           const double* xlower, const double* ylower,
                                           const double* dx, const double* dy,
                                           double rhs[]);

typedef  void (*fc2d_thunderegg_fort_beta_t)(const double* x,
                                            const double* y,
                                            const double* beta,
                                            double grad[]);

typedef double (*fc2d_thunderegg_fort_eval_bc_t)(const int *iface, const double *t, 
                                                const double *x, const double *y);

typedef void (*fc2d_thunderegg_fort_apply_bc_t)(const int* blockno, const  int* mx, const  int* my, 
                                               const  int* mbc, const  int* meqn, 
                                               const double* xlower, const double* ylower,
                                               const double* dx, const double* dy, 
                                               const double *t, 
                                               int intersects_bc[], int mthbc[], 
                                               double rhs[], fc2d_thunderegg_fort_eval_bc_t g_bc, 
                                               int* cons_check, double flux_sum[]);

typedef void (*fc2d_thunderegg_operator_t)(struct fclaw2d_global *glob);

/* -------------------------- Solver and utilities ------------------------------------ */

//void fc2d_thunderegg_solve(struct fclaw2d_global *glob);

/* --------------------------------- Virtual table ------------------------------------ */

struct fc2d_thunderegg_vtable
{
    /* Solver that defines patch operator and calls ThunderEgg solver */
    fc2d_thunderegg_operator_t        patch_operator;  /* 'operator' is a keyword */

	/* Fortran routines */
	fc2d_thunderegg_fort_rhs_t        fort_rhs;	
	fc2d_thunderegg_fort_beta_t       fort_beta;	
    fc2d_thunderegg_fort_apply_bc_t   fort_apply_bc;
    fc2d_thunderegg_fort_eval_bc_t    fort_eval_bc;

	int is_set;
};

/**
 * @brief Initialize the thunderegg solver
 * 
 * fclaw2d_vtables_intialize should be called before this function.
 * 
 * @param global the global context
 */
void fc2d_thunderegg_solver_initialize(struct fclaw2d_global* glob);

/**
 * @brief Get the thunderegg vtable
 * 
 * @param global the global context
 * @return fc2d_thunderegg_vtable_t* the vtable
 */
fc2d_thunderegg_vtable_t* fc2d_thunderegg_vt(struct fclaw2d_global* glob);


/* ----------------------------- User access to solver functions ---------------------- */

void fc2d_thunderegg_setprob(struct fclaw2d_global* glob);


void fc2d_thunderegg_rhs(struct fclaw2d_global* glob,
                        struct fclaw_patch *patch,
                        int blockno,
                        int patchno);


/* -------------------------------- Operator utilities -------------------------------- */

/* Put this here so that user does not have to include fc2d_thunderegg_heat.h */
void fc2d_thunderegg_heat_set_lambda(double lambda);

double fc2d_thunderegg_heat_get_lambda();



#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif /* !FC2D_THUNDEREGG_H */
