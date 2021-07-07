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

#ifndef TRANSPORT_USER_H
#define TRANSPORT_USER_H

#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch_pillow.h>

/* FORTRAN headers are not needed here, but generally needed in the
   advection examples and are included here for convenience */
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_fort.h>

#include <fclaw2d_clawpatch46_fort.h>
#include <fclaw2d_clawpatch5_fort.h>

#include <fc2d_clawpack46.h>  
#include <fc2d_clawpack46_options.h>

#include <fc2d_clawpack5.h>
#include <fc2d_clawpack5_options.h>

#include <fc2d_clawpack46_fort.h>  
#include <fc2d_clawpack5_fort.h>

#include <clawpack46_user_fort.h>  
#include <clawpack5_user_fort.h>

#include "transport_user_fort.h"

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* Only include headers above */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif