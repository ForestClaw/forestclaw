/*
  Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FCLAW_PACKAGE_H
#define FCLAW_PACKAGE_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif
#if defined(__GNUC__) || defined(__clang__)
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif
struct fclaw_global;

/* Opaque pointers */
typedef struct fclaw_package_container fclaw_package_container_t;
typedef struct fclaw_package           fclaw_package_t;

fclaw_package_container_t *fclaw_package_container_new (void);
void fclaw_package_container_destroy (fclaw_package_container_t * pkgs);
int fclaw_package_container_add (fclaw_package_container_t * pkg_container,
                                 void *opt);

/*********************** CODE BELOW STILL USING APP ********************/

/* Create, destroy and add packages */
void fclaw_package_container_new_app (fclaw_app_t *app);
void fclaw_package_container_destroy_app (fclaw_app_t *app);

DEPRECATED int fclaw_package_container_add_pkg(struct fclaw_global* glob,
                                    void* opt);
DEPRECATED void* fclaw_package_get_options(struct fclaw_global *glob, int id);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
