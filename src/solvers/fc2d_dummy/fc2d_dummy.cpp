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

#include "fc2d_dummy.H"

static int dummy_package_id = -1;


/* Patch data is stored in each ClawPatch */
struct patch_data
{
    int value;
};

#if 0
static patch_data*
get_patch_data(ClawPatch *cp)
{
    patch_data *wp =
        (patch_data*) cp->clawpack_patch_data(dummy_package_id);
    return wp;
}
#endif

static void*
patch_data_new()
{
    patch_data* data;
    data = new patch_data;
    return (void*) data;
}

static void
patch_data_delete(void *data)
{
    patch_data *pd = (patch_data*) data;
    FCLAW_ASSERT(pd != NULL);
    delete pd;
}

static const fclaw_package_vtable_t dummy_patch_vtable = {
    patch_data_new,
    patch_data_delete
};


/* -----------------------------------------------------------
   Public interface to routines in this file
   ----------------------------------------------------------- */
int fc2d_dummy_package_register(fclaw_package_container_t *pkg_container,
                                      void* opt)
{
    int id;

    /* Don't register a package more than once */
    FCLAW_ASSERT(dummy_package_id == -1);

    /* Register packages */
    id = fclaw_package_collection_add_pkg(pkg_container,opt,
                                          &dummy_patch_vtable);
    dummy_package_id = id;
    return id;
}
