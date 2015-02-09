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


typedef struct fclaw_package_container fclaw_package_container_t;
typedef struct fclaw_package_data      fclaw_package_data_t;

typedef struct fclaw_package_vtable fclaw_package_vtable_t;
typedef struct fclaw_package        fclaw_package_t;

typedef void* (*fclaw_package_data_new_t)();
typedef void (*fclaw_package_data_delete_t)(void *data);

#define FCLAW_MAX_PACKAGES 20

struct fclaw_package_container
{
    fclaw_package_t *pkgs[FCLAW_MAX_PACKAGES];  /* Make adding packages easy ... */
    int count;
};

struct fclaw_package_vtable
{
    fclaw_package_data_new_t new_patch_data;
    fclaw_package_data_delete_t destroy_patch_data;
};


/* Data attached to a patch that is independent of  what would normally get
   passed in through a parameter call */
struct fclaw_package
{
    void* options;
    fclaw_package_vtable_t vt;
    int id;
};

/* Data associated with each new patch.  A new one of these things will
   get created each time a new ClawPatch is created */
struct fclaw_package_data
{
    void *data[FCLAW_MAX_PACKAGES];
    int count;
};

void* fclaw_package_get_data(fclaw_package_data_t *data_container, int id);


fclaw_package_container_t* fclaw_package_collection_init();

int fclaw_package_collection_add_pkg(fclaw_package_container_t* pkg_container,
                                      void* opt,
                                      const fclaw_package_vtable_t *vtable);

void fclaw_package_collection_destroy(fclaw_package_container_t *pkg_container);


void fclaw_package_patch_data_new(fclaw_package_container_t* pkg_container,
                                  fclaw_package_data_t *data_container);


void fclaw_package_patch_data_destroy(fclaw_package_container_t* pkg_container,
                                      fclaw_package_data_t *data_container);


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
