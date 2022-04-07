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

#ifndef FCLAW2D_VTABLE_H
#define FCLAW2D_VTABLE_H

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw2d_global;
struct fclaw2d_patch;


void fclaw2d_after_regrid(struct fclaw2d_global *glob);

void fclaw2d_after_init(struct fclaw2d_global *glob);

void fclaw_fpe_handling(struct fclaw2d_global *glob);


/* ---------------------------------- Typedefs ---------------------------------------- */  
typedef void (*fclaw2d_vtable_initialize_t)();

typedef void (*fclaw2d_problem_setup_t)(struct fclaw2d_global *glob);

typedef void (*fclaw2d_output_frame_t)(struct fclaw2d_global * glob, int iframe);

typedef void (*fclaw2d_after_regrid_t)(struct fclaw2d_global *glob);

typedef void (*fclaw2d_after_initialization_t)(struct fclaw2d_global *glob);

typedef void (*fclaw_fpe_handling_t)(struct fclaw2d_global *glob);

typedef void (*fclaw_fpe_signal_handler_t)(int sig);

/* ------------------------------------ vtable ---------------------------------------- */  
typedef struct fclaw2d_vtable
{

	fclaw2d_vtable_initialize_t      vtable_init;
	fclaw2d_problem_setup_t          problem_setup;
	fclaw2d_after_initialization_t   after_init;

	/* regridding functions */
	fclaw2d_after_regrid_t           after_regrid;

	/* Output functions */
	fclaw2d_output_frame_t           output_frame;

    fclaw_fpe_handling_t             fpe_enable;
    fclaw_fpe_signal_handler_t       fpe_signal_handler;

	int is_set;


} fclaw2d_vtable_t;


fclaw2d_vtable_t* fclaw2d_vt();

void fclaw2d_vtable_initialize();


#ifdef __cplusplus
}
#endif

#endif
