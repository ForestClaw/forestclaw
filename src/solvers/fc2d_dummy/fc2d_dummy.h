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

#ifndef FC2D_DUMMY_H
#define FC2D_DUMMY_H

#include <fclaw_package.h>
#include <fclaw2d_forestclaw.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

void fc2d_dummy_register(fclaw_app_t *app);

int fc2d_dummy_package_id();

void fc2d_dummy_setup_patch(fclaw_domain_t *domain,
                            fclaw_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx);

void fc2d_dummy_define_data(fclaw_domain_t* domain,fclaw_patch_t* this_patch);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif


#endif
