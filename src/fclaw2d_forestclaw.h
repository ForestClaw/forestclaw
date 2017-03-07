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

#ifndef FCLAW2D_FORESTCLAW_H
#define FCLAW2D_FORESTCLAW_H

/* This file should come first to retrieve configure-time information */
#include <fclaw_config.h>

/* Use as an alternate to GNU feenableexcept */
#ifndef FCLAW_HAVE_FEENABLEEXCEPT
#include <fp_exception_glibc_extension.h>
#endif

#include <fenv.h>
#include <signal.h>

#ifdef FCLAW_HAVE_UNISTD_H
#include <unistd.h>    /* To get process ids */
#endif

/* Basic objects */
#include "forestclaw2d.h"
#include "fclaw2d_base.h"
#include "fclaw2d_domain.h"
#include "fclaw2d_block.h"
/* Don't include fclaw2d_patch.H, since it requires C++ linkage */

/* Basic header files that are probably required by all apps */
#include "fclaw_options.h"
#include "fclaw_register.h"
#include "fclaw_package.h"
#include "fclaw2d_vtable.h"

/* Mapping interface - needed even if 'nomap' is used.  */
#include "p4est_connectivity.h"
#include "fclaw2d_map.h"
#include "fclaw2d_map_query.h"

#include "fclaw_timer.h"

/* These still need to be cleaned up */
#include "fclaw2d_defs.h"

/* Misc */
#include "fclaw_math.h"

#include <fclaw2d_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

void fclaw2d_initialize(fclaw2d_global_t *glob);
void fclaw2d_run(fclaw2d_global_t *glob);
void fclaw2d_finalize(fclaw2d_domain_t **domain);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif


/* Note: either we make this a C .h file, or we remove the extern "C". */

#endif
