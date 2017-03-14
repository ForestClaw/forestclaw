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

struct fclaw2d_global;

/* Opaque pointers */
typedef struct fclaw_package_container fclaw_package_container_t;
typedef struct fclaw_package_data      fclaw_package_data_t;
typedef struct fclaw_package           fclaw_package_t;

/* Needs to be completely defined here */
typedef void* (*fclaw_package_data_new_t)();
typedef void (*fclaw_package_data_delete_t)(void *data);

typedef struct fclaw_package_vtable fclaw_package_vtable_t;

struct fclaw_package_vtable
{
    fclaw_package_data_new_t new_patch_data;
    fclaw_package_data_delete_t destroy_patch_data;
};

fclaw_package_container_t *fclaw_package_container_new (void);
void fclaw_package_container_destroy (fclaw_package_container_t * pkgs);
int fclaw_package_container_add (fclaw_package_container_t * pkg_container,
                                 void *opt,
                                 const fclaw_package_vtable_t *vtable);

/*********************** CODE BELOW STILL USING APP ********************/

/* Create, destroy and add packages */
void fclaw_package_container_new_app (fclaw_app_t *app);
void fclaw_package_container_destroy_app (fclaw_app_t *app);

int fclaw_package_container_add_pkg(fclaw_app_t* app,
                                    void* opt,
                                    const fclaw_package_vtable_t *vtable);

int fclaw_package_container_add_pkg_new(struct fclaw2d_global* glob,
                                        void* opt,
                                        const fclaw_package_vtable_t *vtable);

/* Storage in ClawPatch for data from each package */
void fclaw_package_data_destroy(fclaw_package_data_t* pkg_data);
fclaw_package_data_t* fclaw_package_data_new();
  

/* Create, destroy and add patch data for each package */
void fclaw_package_patch_data_new(fclaw_app_t* app,
                                  fclaw_package_data_t *data_container);


void fclaw_package_patch_data_destroy(fclaw_app_t* app,
                                      fclaw_package_data_t *data_container);

void* fclaw_package_get_data(fclaw_package_data_t *data_container, int id);

void* fclaw_package_get_options(fclaw_app_t* app, int id);
void* fclaw_package_get_options_new(struct fclaw2d_global *glob, int id);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
