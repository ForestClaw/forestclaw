/*
Copyright (c) 2012-2026 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_vtk_vtable.h>
#include <fclaw_physical_bc.h>

static
fclaw_clawpatch_vtk_vtable_entry_t* vtk_vt_entry_new(const char* name,
                                                     fclaw_vtk_entry_type_t type,
                                                     int number_of_components,
                                                     fclaw_vtk_patch_elements_cb_t elements_in_patch,
                                                     fclaw_vtk_patch_cb_t callback)
{
	fclaw_clawpatch_vtk_vtable_entry_t* entry =
		FCLAW_ALLOC(fclaw_clawpatch_vtk_vtable_entry_t, 1);

	entry->name = name;
	entry->type = type;
	entry->number_of_components = number_of_components;
    entry->elements_in_patch = elements_in_patch;
	entry->callback = callback;

	return entry;
}

static
void vtk_vt_entry_init(fclaw_clawpatch_vtk_vtable_entry_t* entry,
                       const char* name,
                       fclaw_vtk_entry_type_t type,
                       int number_of_components,
                       fclaw_vtk_patch_elements_cb_t elements_in_patch,
                       fclaw_vtk_patch_cb_t callback)
{
    entry->name = name;
    entry->type = type;
    entry->number_of_components = number_of_components;
    entry->elements_in_patch = elements_in_patch;
    entry->callback = callback;
}

static
void vtk_add_field_entry(fclaw_clawpatch_vtk_vtable_t* vtk_vt,
						 const char* name,
						 fclaw_vtk_entry_type_t type,
                         int number_of_components,
                         fclaw_vtk_patch_elements_cb_t elements_in_patch,
						 fclaw_vtk_patch_cb_t callback)
{
	fclaw_clawpatch_vtk_vtable_entry_t* entry =
		vtk_vt_entry_new(name, type, number_of_components, elements_in_patch, callback);
	sc_list_append(vtk_vt->field_entries, entry);
}

static
void vtk_add_celldata_entry(fclaw_clawpatch_vtk_vtable_t* vtk_vt,
                            const char* name,
                            fclaw_vtk_entry_type_t type,
                            int number_of_components,
                            fclaw_vtk_patch_elements_cb_t elements_in_patch,
                            fclaw_vtk_patch_cb_t callback)
{
    fclaw_clawpatch_vtk_vtable_entry_t* entry =
        vtk_vt_entry_new(name, type, number_of_components, elements_in_patch, callback);
    sc_list_append(vtk_vt->celldata_entries, entry);
}

static
fclaw_clawpatch_vtk_vtable_t* vtk_vt_new(void)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt =
		(fclaw_clawpatch_vtk_vtable_t*)
		FCLAW_ALLOC_ZERO(fclaw_clawpatch_vtk_vtable_t, 1);

	vtk_vt->field_entries = sc_list_new(NULL);
    vtk_vt->celldata_entries = sc_list_new(NULL);

	return vtk_vt;
}

static
void vtk_vt_destroy(void* vt)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt =
		(fclaw_clawpatch_vtk_vtable_t*) vt;
	fclaw_clawpatch_vtk_vtable_entry_t* entry;

	if (vtk_vt->field_entries != NULL)
	{
		while (vtk_vt->field_entries->elem_count > 0)
		{
			entry = (fclaw_clawpatch_vtk_vtable_entry_t*)
				sc_list_pop(vtk_vt->field_entries);
			FCLAW_FREE(entry);
		}
		sc_list_destroy(vtk_vt->field_entries);
	}
    if (vtk_vt->celldata_entries != NULL)
    {
        while (vtk_vt->celldata_entries->elem_count > 0)
        {
            entry = (fclaw_clawpatch_vtk_vtable_entry_t*)
                sc_list_pop(vtk_vt->celldata_entries);
            FCLAW_FREE(entry);
        }
        sc_list_destroy(vtk_vt->celldata_entries);
    }

	FCLAW_FREE(vt);
}

// default callbacks

// fields

/* Extract physical grid geometry for a patch; works for both 2D and 3D. */
static void
get_patch_geometry (fclaw_global_t * glob, fclaw_patch_t * patch,
                    double *xlower, double *ylower, double *zlower,
                    double *dx, double *dy, double *dz)
{
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        int mx, my, mbc;
        fclaw_clawpatch_2d_grid_data (glob, patch, &mx, &my, &mbc,
                                      xlower, ylower, dx, dy);
        *zlower = 0.0;
        *dz = 0.0;
    }
    else
    {
        int mx, my, mz, mbc;
        fclaw_clawpatch_3d_grid_data (glob, patch, &mx, &my, &mz, &mbc,
                                      xlower, ylower, zlower, dx, dy, dz);
    }
}

static void
write_level_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                int blockno, int patchno,
                fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int32_t *idata = (int32_t *) buffer;
    *idata = (int32_t) patch->level;
}

static void
write_patch_starts_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                       int blockno, int patchno,
                       fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    double xlower, ylower, zlower, dx, dy, dz;
    get_patch_geometry (glob, patch, 
                        &xlower, &ylower, &zlower, &dx, &dy, &dz);
    double *ddata = (double *) buffer;
    *ddata++ = xlower;
    *ddata++ = ylower;
    *ddata++ = zlower;
}

static void
write_patch_spacings_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                         int blockno, int patchno,
                         fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    double xlower, ylower, zlower, dx, dy, dz;
    get_patch_geometry (glob, patch, 
                        &xlower, &ylower, &zlower, &dx, &dy, &dz);
    double *ddata = (double *) buffer;
    *ddata++ = dx;
    *ddata++ = dy;
    *ddata++ = dz;
}

// points

static void
write_2d_coordinate_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
                     int blockno, int patchno,
                     fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    fclaw_map_context_t* cont = fclaw_map_get(glob);

    /* Enumerate point coordinates in the patch */
    double *d = (double *) buffer;
    int jstart,jend,istart,iend;
    if(clawpatch_opt->vtk_ghost_out)
    {
        jstart = -mbc;
        jend = my + mbc;
        istart = -mbc;
        iend = mx + mbc;
    }
    else
    {
        jstart = 0;
        jend = my;
        istart = 0;
        iend = mx;
    }

    for (int j = jstart; j <= jend; ++j)
    {
        const double y = ylower + j * dy;
        for (int i = istart; i <= iend; ++i)
        {
            const double x = xlower + i * dx;
            if (fclaw_opt->manifold)
            {
                double xpp,ypp,zpp;
                FCLAW_MAP_2D_C2M(&cont,&blockno,&x,&y,&xpp,&ypp,&zpp);
                *d++ = xpp;
                *d++ = ypp;
                *d++ = zpp;
            }
            else
            {
                *d++ = x;
                *d++ = y;
                *d++ = 0;
            }
        }
    }
}

static void
write_coordinate_3d_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
                        int blockno, int patchno,
                        fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    fclaw_map_context_t* cont = fclaw_map_get(glob);
    /* Enumerate point coordinates in the patch */
    double *d = (double *) buffer;

    int kstart,kend,jstart,jend,istart,iend;
    if(clawpatch_opt->vtk_ghost_out)
    {
        kstart = -mbc;
        kend = mz + mbc;
        jstart = -mbc;
        jend = my + mbc;
        istart = -mbc;
        iend = mx + mbc;
    }
    else
    {
        kstart = 0;
        kend = mz;
        jstart = 0;
        jend = my;
        istart = 0;
        iend = mx;
    }
    for (int k = kstart; k <= kend; ++k)
    {
        const double z = zlower + k * dz;
        for (int j = jstart; j <= jend; ++j)
        {
            const double y = ylower + j * dy;
            for (int i = istart; i <= iend; ++i)
            {
                const double x = xlower + i * dx;
                if (fclaw_opt->manifold)
                {
                    double xpp,ypp,zpp;
                    FCLAW_MAP_3D_C2M(&cont,&blockno,&x,&y,&z,&xpp,&ypp,&zpp);
                    *d++ = xpp;
                    *d++ = ypp;
                    *d++ = zpp;
                }
                else
                {
                    *d++ = x;
                    *d++ = y;
                    *d++ = z;
                }
            }
        }
    }
}

static void
write_coordinate_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
                    int blockno, int patchno,
                    fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        write_2d_coordinate_cb(glob, patch, blockno, patchno, ctx, buffer);
    }
    else
    {
        write_coordinate_3d_cb(glob, patch, blockno, patchno, ctx, buffer);
    }
}

// cells

static void
write_2d_connectivity_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                          int blockno, int patchno,
                          fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int mx,my,mbc;
    double xlower,ylower,dx,dy;

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                 &xlower,&ylower,&dx,&dy);


    int64_t pbefore;
    if(clawpatch_opt->vtk_ghost_out)
    {
        pbefore = (mx + 1 + 2*mbc) * (my + 1 + 2*mbc) *
            (glob->domain->global_num_patches_before +
             glob->domain->blocks[blockno].num_patches_before + patchno);
    }
    else
    {
        pbefore = (mx + 1) * (my + 1) *
            (glob->domain->global_num_patches_before +
             glob->domain->blocks[blockno].num_patches_before + patchno);
    }

    int iend,jend, ystride;
    if(clawpatch_opt->vtk_ghost_out)
    {
        iend = mx + 2*mbc;
        jend = my + 2*mbc;
        ystride = mx + 1 + 2*mbc;
    }
    else
    {
        iend = mx;
        jend = my;
        ystride = mx + 1;
    }
    if (ctx->fits32)
    {
        int32_t *idata = (int32_t *) buffer;
        int32_t l;

        for (int j = 0; j < jend; ++j)
        {
            for (int i = 0; i < iend; ++i)
            {
                l = (int32_t) pbefore + i + j * ystride;
                *idata++ = l;
                *idata++ = l + 1;
                *idata++ = l + ystride + 1;
                *idata++ = l + ystride;
            }
        }
    }
    else
    {
        int64_t *idata = (int64_t *) buffer;
        int64_t l;
        for (int j = 0; j < jend; ++j)
        {
            for (int i = 0; i < iend; ++i)
            {
                l = pbefore + i + j * ystride;
                *idata++ = l;
                *idata++ = l + 1;
                *idata++ = l + ystride + 1;
                *idata++ = l + ystride;
            }
        }
    }
}

static void
write_3d_connectivity_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                          int blockno, int patchno,
                          fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    fclaw_clawpatch_3d_grid_data (glob, patch, &mx, &my, &mz, &mbc,
                                  &xlower, &ylower, &zlower,
                                  &dx, &dy, &dz);

    int64_t pbefore;
    if(clawpatch_opt->vtk_ghost_out)
    {
        pbefore = (mx + 1 + 2*mbc) * (my + 1 + 2*mbc) * (mz + 1 + 2*mbc) *
            (glob->domain->global_num_patches_before +
             glob->domain->blocks[blockno].num_patches_before + patchno);
    }
    else
    {
        pbefore = (mx + 1) * (my + 1) * (mz + 1) *
            (glob->domain->global_num_patches_before +
             glob->domain->blocks[blockno].num_patches_before + patchno);
    }

    int iend, jend, kend, ystride, zstride;
    if(clawpatch_opt->vtk_ghost_out)
    {
        iend = mx + 2*mbc;
        jend = my + 2*mbc;
        kend = mz + 2*mbc;
        ystride = mx + 1 + 2*mbc;
        zstride = (mx + 1 + 2*mbc) * (my + 1 + 2*mbc);
    }
    else
    {
        iend = mx;
        jend = my;
        kend = mz;
        ystride = mx + 1;
        zstride = (mx + 1) * (my + 1);
    }

    if (ctx->fits32)
    {
        int32_t *idata = (int32_t *) buffer;
        int32_t l;

        for (int k = 0; k < kend; ++k)
        {
            for (int j = 0; j < jend; ++j)
            {
                for (int i = 0; i < iend; ++i)
                {
                    l = (int32_t) pbefore + i + j * ystride
                         + k * zstride;
                    *idata++ = l;
                    *idata++ = l + 1;
                    *idata++ = l + ystride + 1;
                    *idata++ = l + ystride;
                    *idata++ = l + zstride;
                    *idata++ = l + zstride + 1;
                    *idata++ = l + zstride + ystride + 1;
                    *idata++ = l + zstride + ystride;
                }
            }
        }
    }
    else
    {
        int64_t *idata = (int64_t *) buffer;
        int64_t l;

        for (int k = 0; k < kend; ++k)
        {
            for (int j = 0; j < jend; ++j)
            {
                for (int i = 0; i < iend; ++i)
                {
                    l = pbefore + i + j * ystride
                         + k * zstride;
                    *idata++ = l;
                    *idata++ = l + 1;
                    *idata++ = l + ystride + 1;
                    *idata++ = l + ystride;
                    *idata++ = l + zstride;
                    *idata++ = l + zstride + 1;
                    *idata++ = l + zstride + ystride + 1;
                    *idata++ = l + zstride + ystride;
                }
            }
        }
    }
}

static void
write_connectivity_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                       int blockno, int patchno,
                       fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        write_2d_connectivity_cb(glob, patch, blockno, patchno, ctx, buffer);
    }
    else
    {
        write_3d_connectivity_cb(glob, patch, blockno, patchno, ctx, buffer);
    }
}

static int
get_cells_per_patch(fclaw_global_t *glob, fclaw_patch_t *patch)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;

    int num_cells = 0;
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        fclaw_clawpatch_2d_grid_data(glob, patch, &mx, &my, &mbc,
                                     &xlower, &ylower, &dx, &dy);

        if(clawpatch_opt->vtk_ghost_out)
        {
            num_cells = (mx + 2*mbc) * (my + 2*mbc);
        }
        else
        {
            num_cells = mx * my;
        }
    }
    else
    {
        fclaw_clawpatch_3d_grid_data(glob, patch, &mx, &my, &mz, &mbc,
                                     &xlower, &ylower, &zlower,
                                     &dx, &dy, &dz);
        if(clawpatch_opt->vtk_ghost_out)
        {
            num_cells = (mx + 2*mbc) * (my + 2*mbc) * (mz + 2*mbc);
        }
        else
        {
            num_cells = mx * my * mz;
        }
    }

    return num_cells;
}

static void
write_offsets_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                  int blockno, int patchno,
                  fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int num_points_per_cell;
    if(fclaw_clawpatch_dim(patch) == 2)
    {
        num_points_per_cell = 4;
    }
    else
    {
        num_points_per_cell = 8;
    }

    const int cells_per_patch = get_cells_per_patch(glob, patch);
    const int64_t cbefore = cells_per_patch *
        (glob->domain->global_num_patches_before +
         glob->domain->blocks[blockno].num_patches_before + patchno);

    if (ctx->fits32)
    {
        int32_t *idata = (int32_t *) buffer;
        if(ctx->offsets_include_zero && glob->mpirank == 0 && patchno == 0 && blockno == 0)
        {
            *idata++ = 0;
        }
        int32_t k = num_points_per_cell * (int32_t) (cbefore + 1);
        for (int c = 0; c < cells_per_patch; k += num_points_per_cell, ++c)
        {
            *idata++ = k;
        }
    }
    else
    {
        int64_t *idata = (int64_t *) buffer;
        if(ctx->offsets_include_zero && glob->mpirank == 0 && patchno == 0 && blockno == 0)
        {
            *idata++ = 0;
        }
        int64_t k = num_points_per_cell * (cbefore + 1);
        for (int c = 0; c < cells_per_patch; k += num_points_per_cell, ++c)
        {
            *idata++ = k;
        }
    }
}

static void
write_types_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                int blockno, int patchno,
                fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    char *cdata = buffer;
    const int cells_per_patch = get_cells_per_patch(glob, patch);
    if(fclaw_clawpatch_dim(patch) == 2)
    {
        for (int c = 0; c < cells_per_patch; ++c)
        {
            *cdata++ = 9;
        }
    }
    else 
    {
        for (int c = 0; c < cells_per_patch; ++c)
        {
            *cdata++ = 12;
        }
    }
}

// celldata 

static void
write_mpirank_cb (fclaw_global_t* glob,
                  fclaw_patch_t* patch,
                  int blockno, int patchno,
                  fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int *idata = (int *) buffer;
    const int cells_per_patch = get_cells_per_patch(glob, patch);
    for (int c = 0; c < cells_per_patch; ++c)
    {
        *idata++ = glob->mpirank;
    }
}

static void
write_blockno_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                  int blockno, int patchno,
                  fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int *idata = (int *) buffer;
    const int cells_per_patch = get_cells_per_patch(glob, patch);
    for (int c = 0; c < cells_per_patch; ++c)
    {
        *idata++ = blockno;
    }
}

static void
write_patchno_cb (fclaw_global_t* glob, fclaw_patch_t* patch,
                  int blockno, int patchno,
                  fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int *idata = (int32_t *) buffer;
    const int cells_per_patch = get_cells_per_patch(glob, patch);
    for (int c = 0; c < cells_per_patch; ++c)
    {
        *idata++ = patchno;
    }
}

/* typdef for data acess functions q,rhs,etc */
typedef void (*patch_data_access_t)(struct fclaw_global *glob, struct fclaw_patch *patch, double **q, int *meqn);

static void
pack_data (fclaw_global_t * glob,
           fclaw_patch_t * patch,
           int blockno, int patchno,
           patch_data_access_t access,
           char *a)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int meqn;
    double *q;
    access(glob, patch, &q, &meqn);

    int patch_dim = fclaw_clawpatch_dim(patch);
    if(patch_dim == 2)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;

        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        const int xlane = mx + 2 * mbc;
        const int ylane = my + 2 * mbc;

        // Enumerate equation data in the patch
        float *f = (float *) a;
        int istart,iend,jstart,jend;
        if(clawpatch_opt->vtk_ghost_out)
        {
            jstart = -mbc;
            jend = my + mbc;
            istart = -mbc;
            iend = mx + mbc;
        }
        else
        {
            jstart = 0;
            jend = my;
            istart = 0;
            iend = mx;
        }
        for (int j = jstart; j < jend; ++j)
        {
            for (int i = istart; i < iend; ++i)
            {
                for (int k = 0; k < meqn; ++k)
                {
                    /* For Clawpack 5.0 layout */
                    //*f++ = (float) q[((j+mbc)*xlane + (i+mbc))*meqn + k];

                    /* For Clawpack 4.x layout */
                    *f++ = (float) q[(k * ylane + j + mbc) * xlane + i + mbc];
                }
            }
        }

    }
    else
    {
        int mx,my,mz,mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                   &xlower,&ylower,&zlower, &dx,&dy, &dz);

        const int xlane = mx + 2 * mbc;
        const int ylane = my + 2 * mbc;
        const int zlane = mz + 2 * mbc;

        // Enumerate equation data in the patch
        float *f = (float *) a;
        int istart,iend,jstart,jend,kstart,kend;
        if(clawpatch_opt->vtk_ghost_out)
        {
            kstart = -mbc;
            kend = mz + mbc;
            jstart = -mbc;
            jend = my + mbc;
            istart = -mbc;
            iend = mx + mbc;
        }
        else
        {
            kstart = 0;
            kend = mz;
            jstart = 0;
            jend = my;
            istart = 0; 
            iend = mx;
        }
        for (int k = kstart; k < kend; ++k)
        {
            for (int j = jstart; j < jend; ++j)
            {
                for (int i = istart; i < iend; ++i)
                {
                    for (int eqn = 0; eqn < meqn; ++eqn)
                    {
                        /* For Clawpack 5.0 layout */
                        //*f++ = (float) q[((j+mbc)*xlane + (i+mbc))*meqn + k];

                        /* For Clawpack 4.x layout */
                        *f++ = (float) q[eqn * zlane * ylane * xlane + (k + mbc) * ylane * xlane + (j + mbc) * xlane + i + mbc];
                    }
                }
            }
        }
    }
}

static void
write_value_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
                int blockno, int patchno,
                fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_soln_data, buffer);
}

static void
write_aux_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
              int blockno, int patchno,
              fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int maux;
    double *q;
    fclaw_clawpatch_aux_data(glob, patch, &q, &maux);

    int patch_dim = fclaw_clawpatch_dim(patch);
    if(patch_dim == 2)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;

        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        const int xlane = mx + 2 * mbc;
        const int ylane = my + 2 * mbc;

        // Enumerate equation data in the patch
        float *f = (float *) buffer;
        int istart,iend,jstart,jend;
        if(clawpatch_opt->vtk_ghost_out)
        {
            jstart = -mbc;
            jend = my + mbc;
            istart = -mbc;
            iend = mx + mbc;
        }
        else
        {
            jstart = 0;
            jend = my;
            istart = 0;
            iend = mx;
        }
        for (int j = jstart; j < jend; ++j)
        {
            for (int i = istart; i < iend; ++i)
            {
                for (int k = 0; k < maux; ++k)
                {
                    if(clawpatch_opt->vtk_aux_out[k] > 0)
                    {
                        /* For Clawpack 5.0 layout */
                        //*f++ = (float) q[((j+mbc)*xlane + (i+mbc))*meqn + k];

                        /* For Clawpack 4.x layout */
                        int aux_out = clawpatch_opt->vtk_aux_out[k] - 1;
                        *f++ = (float) q[(aux_out * ylane + j + mbc) * xlane + i + mbc];
                    }
                }
            }
        }

    }
    else
    {
        int mx,my,mz,mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                   &xlower,&ylower,&zlower, &dx,&dy, &dz);

        const int xlane = mx + 2 * mbc;
        const int ylane = my + 2 * mbc;
        const int zlane = mz + 2 * mbc;

        // Enumerate equation data in the patch
        float *f = (float *) buffer;
        int istart,iend,jstart,jend,kstart,kend;
        if(clawpatch_opt->vtk_ghost_out)
        {
            kstart = -mbc;
            kend = mz + mbc;
            jstart = -mbc;
            jend = my + mbc;
            istart = -mbc;
            iend = mx + mbc;
        }
        else
        {
            kstart = 0;
            kend = mz;
            jstart = 0;
            jend = my;
            istart = 0;
            iend = mx;
        }
        for (int k = kstart; k < kend; ++k)
        {
            for (int j = jstart; j < jend; ++j)
            {
                for (int i = istart; i < iend; ++i)
                {
                    for (int eqn = 0; eqn < maux; ++eqn)
                    {
                        if(clawpatch_opt->vtk_aux_out[eqn])
                        {
                            /* For Clawpack 5.0 layout */
                            //*f++ = (float) q[((j+mbc)*xlane + (i+mbc))*meqn + k];

                            /* For Clawpack 4.x layout */
                            int aux_out = clawpatch_opt->vtk_aux_out[eqn] - 1;
                            *f++ = (float) q[aux_out * zlane * ylane * xlane + (k + mbc) * ylane * xlane + (j + mbc) * xlane + i + mbc];
                        }
                    }
                }
            }
        }
    }

}

static void
write_rhs_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
              int blockno, int patchno,
              fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_rhs_data, buffer);
}

static void
write_soln_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
               int blockno, int patchno,
               fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_elliptic_soln_data, buffer);
}

static void
write_error_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
                int blockno, int patchno,
                fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_elliptic_error_data, buffer);
}

static void
write_ghost_cb (fclaw_global_t * glob, fclaw_patch_t * patch,
                int blockno, int patchno,
                fclaw_vtk_cb_context_t* ctx, char* buffer)
{
    int patch_dim = fclaw_clawpatch_dim(patch);
    if(patch_dim == 2)
    {
        int mx,my,mbc;
        double xlower,ylower,dx,dy;

        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);

        // Enumerate equation data in the patch
        for (int j = -mbc; j < my + mbc; ++j)
        {
            for (int i = -mbc; i < mx + mbc; ++i)
            {
                if(i < 0 || i >= mx || j < 0 || j >= my)
                {
                    *buffer ++ = 32;
                }
                else
                {
                    *buffer ++ = 0;
                }
            }
        }

    }
    else
    {
        int mx,my,mz,mbc;
        double xlower,ylower,zlower,dx,dy,dz;
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                   &xlower,&ylower,&zlower, &dx,&dy, &dz);

        //get physcial boundaries
        int intersects_bc[6] = {1,1,1,1,1,1};
        fclaw_physical_get_bc(glob, blockno, patchno, intersects_bc);
        for(int k = -mbc; k < mz + mbc; ++k)
        {
            for (int j = -mbc; j < my + mbc; ++j)
            {
                for (int i = -mbc; i < mx + mbc; ++i)
                {
                    if(i < 0 || i >= mx || j < 0 || j >= my || k < 0 || k >= mz)
                    {
                        if((i < 0 && intersects_bc[0]) || (i >= mx && intersects_bc[1]) ||
                           (j < 0 && intersects_bc[2]) || (j >= my && intersects_bc[3]) ||
                           (k < 0 && intersects_bc[4]) || (k >= mz && intersects_bc[5]))
                        {
                            *buffer ++ = 32;
                        }
                        else
                        {
                            *buffer ++ = 32;
                        }
                    }
                    else
                    {
                        *buffer ++ = 0;
                    }
                }
            }
        }
    }
}

// size callbacks

static size_t
one_element_per_patch (fclaw_global_t* glob, fclaw_patch_t* patch, 
                       int blockno, int patchno, 
                       fclaw_vtk_cb_context_t* ctx)
{
    return 1;
}

static size_t
points_in_patch (fclaw_global_t* glob, fclaw_patch_t* patch, 
                 int blockno, int patchno,
                 fclaw_vtk_cb_context_t* ctx)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;
    size_t num_points = 0;
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        fclaw_clawpatch_2d_grid_data(glob, patch, &mx, &my, &mbc, 
                                     &xlower, &ylower, &dx, &dy);
        if(clawpatch_opt->vtk_ghost_out)
        {
            num_points = (size_t) (mx + 1 + 2*mbc) * (my + 1 + 2*mbc);
        }
        else
        {
            num_points = (size_t) (mx + 1) * (my + 1);
        }
    }
    else
    {
        fclaw_clawpatch_3d_grid_data(glob, patch, &mx, &my, &mz, &mbc, 
                                     &xlower, &ylower, &zlower, &dx, &dy, &dz);
        if(clawpatch_opt->vtk_ghost_out)
        {
            num_points = (size_t) (mx + 1 + 2*mbc) * (my + 1 + 2*mbc) * (mz + 1 + 2*mbc);
        }
        else
        {
            num_points = (size_t) (mx + 1) * (my + 1) * (mz + 1);
        }
    }
    return num_points;
}

static size_t
cells_in_patch (fclaw_global_t* glob, fclaw_patch_t* patch, 
                int blockno, int patchno,
                fclaw_vtk_cb_context_t* ctx)
{
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;

    int num_cells = 0;
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        fclaw_clawpatch_2d_grid_data(glob, patch, &mx, &my, &mbc, 
                                     &xlower, &ylower, &dx, &dy);
        if(clawpatch_opt->vtk_ghost_out)
        {
            num_cells = (size_t) (mx + 2*mbc) * (my + 2*mbc);
        }
        else
        {
            num_cells = (size_t) mx * my;
        }
    }
    else
    {
        fclaw_clawpatch_3d_grid_data(glob, patch, &mx, &my, &mz, &mbc, 
                                     &xlower, &ylower, &zlower, &dx, &dy, &dz);
        if(clawpatch_opt->vtk_ghost_out)
        {
            num_cells = (size_t) (mx + 2*mbc) * (my + 2*mbc) * (mz + 2*mbc);
        }
        else
        {
            num_cells = (size_t) mx * my * mz;
        }
    }
    return num_cells;
}

static size_t
connectivity_in_patch (fclaw_global_t* glob, fclaw_patch_t* patch, 
                       int blockno, int patchno,
                       fclaw_vtk_cb_context_t* ctx)
{
    int num_cells = cells_in_patch(glob, patch, blockno, patchno, ctx);
    int num_points_per_cell;
    if(fclaw_clawpatch_dim(patch) == 2)
    {
        num_points_per_cell = 4;
    }
    else
    {
        num_points_per_cell = 8;
    }
    return (size_t) num_cells * num_points_per_cell;
}

static size_t
offsets_in_patch (fclaw_global_t* glob, fclaw_patch_t* patch, 
                  int blockno, int patchno,
                  fclaw_vtk_cb_context_t* ctx)
{
    size_t retval = cells_in_patch(glob, patch, blockno, patchno, ctx);
    if(ctx->offsets_include_zero && glob->mpirank == 0 && patchno == 0 && blockno == 0)
    {
        retval += 1;
    }
    return retval;
}

// Public interface to add entries to the vtable
fclaw_clawpatch_vtk_vtable_t*
fclaw_clawpatch_vtk_vtable(struct fclaw_global* glob)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt =
		(fclaw_clawpatch_vtk_vtable_t*)
		fclaw_global_get_vtable((fclaw_global_t*) glob,
								"fclaw_clawpatch_vtk_vtable");
	FCLAW_ASSERT(vtk_vt != NULL);
	return vtk_vt;
}

void fclaw_clawpatch_vtk_vtable_initialize(struct fclaw_global* glob)
{
	fclaw_clawpatch_vtk_vtable_t* vtk_vt = vtk_vt_new();
    fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    // field entries
    vtk_add_field_entry(vtk_vt, "levels", FCLAW_VTK_INT32, 
                        1,
                        &one_element_per_patch,
                        &write_level_cb);
    vtk_add_field_entry(vtk_vt, "patch_starts", FCLAW_VTK_FLOAT64, 
                        3,
                        &one_element_per_patch,
                        &write_patch_starts_cb);
    vtk_add_field_entry(vtk_vt, "patch_spacings", FCLAW_VTK_FLOAT64, 
                        3,
                        &one_element_per_patch,
                        &write_patch_spacings_cb);

    // primary point entry
    vtk_vt_entry_init(&vtk_vt->position_entry, "Position", FCLAW_VTK_FLOAT64,
                      3,
                      &points_in_patch,
                      &write_coordinate_cb);

    // primary cell entries
    vtk_vt_entry_init(&vtk_vt->connectivity_entry, "connectivity", FCLAW_VTK_INT32_OR_64,
                      1,
                      &connectivity_in_patch,
                      &write_connectivity_cb);
    vtk_vt_entry_init(&vtk_vt->offsets_entry, "offsets", FCLAW_VTK_INT32_OR_64,
                      1,
                      &offsets_in_patch,
                      &write_offsets_cb);
    vtk_vt_entry_init(&vtk_vt->types_entry, "types", FCLAW_VTK_UINT8,
                      1,
                      &cells_in_patch,
                      &write_types_cb);

    int meqn = clawpatch_opt->meqn;
    // celldata entries
    vtk_add_celldata_entry(vtk_vt, "mpirank", FCLAW_VTK_INT32,
                           1,
                           &cells_in_patch,
                           &write_mpirank_cb);
    vtk_add_celldata_entry(vtk_vt, "blockno", FCLAW_VTK_INT32,
                           1,
                           &cells_in_patch,
                           &write_blockno_cb);
    vtk_add_celldata_entry(vtk_vt, "patchno", FCLAW_VTK_INT32,
                           1,
                            &cells_in_patch,
                           &write_patchno_cb);

    vtk_add_celldata_entry(vtk_vt, "meqn", FCLAW_VTK_FLOAT32,
                           meqn,
                           &cells_in_patch,
                           &write_value_cb);
    
    int num_aux_fields = 0;
    for (int i = 0; i < clawpatch_opt->maux; ++i)
    {
        if (clawpatch_opt->vtk_aux_out[i] > 0)
        {
            num_aux_fields++;
        }
    }

    if (num_aux_fields > 0)
    {
        vtk_add_celldata_entry(vtk_vt, "aux", FCLAW_VTK_FLOAT32,
                               num_aux_fields,
                               &cells_in_patch,
                               &write_aux_cb);
    }

    int rhs_fields = clawpatch_opt->rhs_fields;
    if (rhs_fields > 0)
    {
        vtk_add_celldata_entry(vtk_vt, "rhs", FCLAW_VTK_FLOAT32,
                               rhs_fields,
                               &cells_in_patch,
                               &write_rhs_cb);
        vtk_add_celldata_entry(vtk_vt, "soln", FCLAW_VTK_FLOAT32,
                               rhs_fields,
                               &cells_in_patch,
                               &write_soln_cb);
        vtk_add_celldata_entry(vtk_vt, "error", FCLAW_VTK_FLOAT32,
                               rhs_fields,
                               &cells_in_patch,
                               &write_error_cb);
    }

    if(clawpatch_opt->vtk_ghost_out)
    {
        vtk_add_celldata_entry(vtk_vt, "vtkGhostType", FCLAW_VTK_UINT8,
                               1,
                               &cells_in_patch,
                               &write_ghost_cb);
    }

	fclaw_global_vtable_store((fclaw_global_t*) glob,
							  "fclaw_clawpatch_vtk_vtable",
							  vtk_vt,
							  vtk_vt_destroy);
}


