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

#include <fclaw_package.h>


#define FCLAW_MAX_PACKAGES 20


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

struct fclaw_package_container
{
    fclaw_package_t *pkgs[FCLAW_MAX_PACKAGES];  /* Make adding packages easy ... */
    int count;
    int max_packages;
};


fclaw_package_data_t* fclaw_package_data_new()
{
    fclaw_package_data_t* pkg_data = FCLAW_ALLOC(fclaw_package_data_t,1);
    return pkg_data;
}

void fclaw_package_data_destroy(fclaw_package_data_t* pkg_data)
{
    FCLAW_ASSERT(pkg_data != NULL);
    FCLAW_FREE(pkg_data);
    pkg_data = NULL;
}

fclaw_package_container_t *
fclaw_package_container_new (void)
{
    int i;
    fclaw_package_container_t* pkg_container;

    pkg_container = FCLAW_ALLOC(fclaw_package_container_t,1);
    for(i = 0; i < FCLAW_MAX_PACKAGES; i++)
    {
        pkg_container->pkgs[i] = NULL;
    }
    pkg_container->count = 0;
    pkg_container->max_packages = FCLAW_MAX_PACKAGES;

    return pkg_container;
};

void fclaw_package_container_new_app (fclaw_app_t* app)
{
    fclaw_app_set_attribute(app,"packages",
                            fclaw_package_container_new ());
};

void
fclaw_package_container_destroy (fclaw_package_container_t * pkg_container)
{
    int i;
    fclaw_package_t *pkg;

    FCLAW_ASSERT (pkg_container != NULL);

    for(i = 0; i < pkg_container->count; i++)
    {
        pkg = pkg_container->pkgs[i];
        FCLAW_ASSERT(pkg != NULL);
        FCLAW_FREE(pkg);
    }
    FCLAW_FREE(pkg_container);
}

void
fclaw_package_container_destroy_app (fclaw_app_t *app)
{
    fclaw_package_container_t * pkg_container = (fclaw_package_container_t *)
      fclaw_app_get_attribute (app, "packages", NULL);

    fclaw_package_container_destroy (pkg_container);
}

int
fclaw_package_container_add (fclaw_package_container_t * pkg_container,
                             void *opt, const fclaw_package_vtable_t *vtable)
{
    int id;
    fclaw_package_t *new_pkg;

    FCLAW_ASSERT (pkg_container != NULL);

    FCLAW_ASSERT(pkg_container->count < FCLAW_MAX_PACKAGES);
    id = pkg_container->count++;
    new_pkg = FCLAW_ALLOC(fclaw_package_t,1);
    new_pkg->id = id;
    new_pkg->vt = *vtable;
    new_pkg->options = opt;
    pkg_container->pkgs[id] = new_pkg;
    return id;
}

int
    fclaw_package_container_add_pkg(fclaw_app_t* app,
                                    void* opt,
                                    const fclaw_package_vtable_t *vtable)
{
    fclaw_package_container_t *pkg_container;

    pkg_container = (fclaw_package_container_t *)
      fclaw_app_get_attribute (app, "packages", NULL);

    return fclaw_package_container_add (pkg_container, opt, vtable);
}

void
    fclaw_package_patch_data_new(fclaw_app_t* app,
                                 fclaw_package_data_t *data_container)
{
    int i;
    fclaw_package_t *pkg;
    fclaw_package_container_t* pkg_container;

    pkg_container = (fclaw_package_container_t *)
      fclaw_app_get_attribute(app,"packages",NULL);

    if (pkg_container != NULL)
    {

        for(i = 0; i < pkg_container->count; i++)
        {
            pkg = pkg_container->pkgs[i];
            FCLAW_ASSERT(pkg != NULL);
            data_container->data[i] = (void*) pkg->vt.new_patch_data();
        }
        data_container->count = pkg_container->count;
    }
}

void
    fclaw_package_patch_data_destroy(fclaw_app_t *app,
                                     fclaw_package_data_t *data_container)
{
    int i;
    fclaw_package_t *pkg;
    fclaw_package_container_t* pkg_container;

    pkg_container = (fclaw_package_container_t *)
      fclaw_app_get_attribute(app, "packages", NULL);

    if (pkg_container != NULL)
    {
        for(i = 0; i < pkg_container->count; i++)
        {
            pkg = pkg_container->pkgs[i];
            FCLAW_ASSERT(pkg != NULL);
            pkg->vt.destroy_patch_data(data_container->data[i]);
        }
    }
}

void* fclaw_package_get_data(fclaw_package_data_t *data_container, int id)
{
    FCLAW_ASSERT(0 <= id  && id < data_container->count);
    return data_container->data[id];
}

void* fclaw_package_get_options(fclaw_app_t* app, int id)
{
    fclaw_package_t *pkg;
    fclaw_package_container_t* pkg_container;

    pkg_container = (fclaw_package_container_t *)
      fclaw_app_get_attribute (app, "packages", NULL);
    FCLAW_ASSERT (pkg_container != NULL);
    FCLAW_ASSERT (0 <= id && id < pkg_container->count);
    pkg = pkg_container->pkgs[id];
    return pkg->options;
}
