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

#ifndef FCLAW3D_INCLUDE_ALL_H
#define FCLAW3D_INCLUDE_ALL_H

#include <fclaw_base.h>
#include <fclaw_domain.h>
#include <fclaw_options.h>
#include <fclaw_global.h>

#include <fclaw3d_defs.h>

#ifdef P8HACK

#include <fclaw_patch.h>

#include <fclaw3d_forestclaw.h>
#include <fclaw3d_diagnostics.h>

#endif /* P8HACK */

#include <fclaw3d_convenience.h>

#ifdef P8HACK

#include <fclaw3d_vtable.h>

#include <fclaw3d_physical_bc.h>

#include <fclaw_gauges.h>

#include <fclaw3d_map.h>
#include <fclaw3d_map_brick.h>
#include <fclaw3d_map_query.h>

#include <fclaw_math.h>

#endif /* P8HACK */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* This file is for convenience only and should only be used in 
	applications ... */

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !FCLAW3D_INCLUDE_ALL_H */
