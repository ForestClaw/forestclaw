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

#include <fclaw_clawpatch_output_hdf5.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fclaw_global.h>

#include <fclaw_options.h>
#include <fclaw2d_map.h>
#include <hdf5_hl.h>
#include <hdf5.h>

typedef struct fclaw2d_hdf5_state
{
    int dim;
    int patch_children;
    int mx, my, mz;
    int meqn;
    int points_per_patch, cells_per_patch;
    int intsize, ndsize;
    int fits32;
    char filename[BUFSIZ];
    int64_t global_num_points, global_num_cells;
    int64_t global_num_connectivity;
    int64_t offset_position, psize_position;
    int64_t offset_connectivity, psize_connectivity;
    int64_t offset_offsets, psize_offsets;
    int64_t offset_types, psize_types;
    int64_t offset_mpirank, psize_mpirank;
    int64_t offset_blockno, psize_blockno;
    int64_t offset_patchno, psize_patchno;
    int64_t offset_meqn, psize_meqn;
    int64_t offset_end;
    const char *inttype;
    fclaw_hdf5_patch_data_t coordinate_cb;
    fclaw_hdf5_patch_data_t value_cb;
    FILE *file;
#ifdef P4EST_ENABLE_MPIIO
    MPI_File mpifile;
    MPI_Offset mpibegin;
#endif
    char *buf;
}
fclaw2d_hdf5_state_t;

/**
 * @brief Write the buffer to file
 *
 * @param s the hdf5 state
 * @param psize_field the size of the buffer
 */
static void
write_buffer (fclaw2d_hdf5_state_t * s, int64_t psize_field)
{
#ifndef P4EST_ENABLE_MPIIO
    size_t retvalz;

    retvalz = fwrite (s->buf, psize_field, 1, s->file);
    SC_CHECK_ABORT (retvalz == 1, "VTK file write failed");
#else
    int mpiret;
    MPI_Status mpistatus;

    mpiret = MPI_File_write (s->mpifile, s->buf, psize_field, MPI_BYTE,
                             &mpistatus);
    SC_CHECK_MPI (mpiret);
#endif
}

static void
write_position_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                   int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;

    s->coordinate_cb (g->glob, patch, blockno, patchno, s->buf);
    write_buffer (s, s->psize_position);
}

static void
write_2d_connectivity_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                       int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int i, j;
    const int64_t pbefore = s->points_per_patch *
        (domain->global_num_patches_before +
         domain->blocks[blockno].num_patches_before + patchno);

    if (s->fits32)
    {
        int32_t *idata = (int32_t *) s->buf;
        int32_t l;

        for (j = 0; j < s->my; ++j)
        {
            for (i = 0; i < s->mx; ++i)
            {
                l = (int32_t) pbefore + i + j * (s->mx + 1);
                *idata++ = l;
                *idata++ = l + 1;
                *idata++ = l + (s->mx + 2);
                *idata++ = l + (s->mx + 1);
            }
        }
    }
    else
    {
        int64_t *idata = (int64_t *) s->buf;
        int64_t l;
        for (j = 0; j < s->my; ++j)
        {
            for (i = 0; i < s->mx; ++i)
            {
                l = pbefore + i + j * (s->mx + 1);
                *idata++ = l;
                *idata++ = l + 1;
                *idata++ = l + (s->mx + 2);
                *idata++ = l + (s->mx + 1);
            }
        }
    }
    write_buffer (s, s->psize_connectivity);
}
static void
write_3d_connectivity_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                          int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int i, j, k;
    const int64_t pbefore = s->points_per_patch *
        (domain->global_num_patches_before +
         domain->blocks[blockno].num_patches_before + patchno);

    if (s->fits32)
    {
        int32_t *idata = (int32_t *) s->buf;
        int32_t l;

        for (k = 0; k < s->mz; ++k)
        {
            for (j = 0; j < s->my; ++j)
            {
                for (i = 0; i < s->mx; ++i)
                {
                    l = (int32_t) pbefore + i + j * (s->mx + 1)
                         + k * (s->my + 1) * (s->mx + 1);
                    *idata++ = l;
                    *idata++ = l + 1;
                    *idata++ = l + (s->mx + 2);
                    *idata++ = l + (s->mx + 1);
                    *idata++ = l + (s->mx + 1) * (s->my + 1);
                    *idata++ = l + (s->mx + 1) * (s->my + 1) + 1;
                    *idata++ = l + (s->mx + 1) * (s->my + 1) + (s->mx + 2);
                    *idata++ = l + (s->mx + 1) * (s->my + 1) + (s->mx + 1);
                }
            }
        }
    }
    else
    {
        int64_t *idata = (int64_t *) s->buf;
        int64_t l;

        for (k = 0; k < s->mz; ++k)
        {
            for (j = 0; j < s->my; ++j)
            {
                for (i = 0; i < s->mx; ++i)
                {
                    l = pbefore + i + j * (s->mx + 1)
                         + k * (s->my + 1) * (s->mx + 1);
                    *idata++ = l;
                    *idata++ = l + 1;
                    *idata++ = l + (s->mx + 2);
                    *idata++ = l + (s->mx + 1);
                    *idata++ = l + (s->mx + 1) * (s->my + 1);
                    *idata++ = l + (s->mx + 1) * (s->my + 1) + 1;
                    *idata++ = l + (s->mx + 1) * (s->my + 1) + (s->mx + 2);
                    *idata++ = l + (s->mx + 1) * (s->my + 1) + (s->mx + 1);
                }
            }
        }
    }
    write_buffer (s, s->psize_connectivity);
}

static void
write_offsets_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                  int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int c;
    const int64_t cbefore = s->cells_per_patch *
        (domain->global_num_patches_before +
         domain->blocks[blockno].num_patches_before + patchno);

    if (s->fits32)
    {
        int32_t *idata = (int32_t *) s->buf;
        int32_t k = s->patch_children * (int32_t) (cbefore + 1);
        for (c = 0; c < s->cells_per_patch; k += s->patch_children, ++c)
        {
            *idata++ = k;
        }
    }
    else
    {
        int64_t *idata = (int64_t *) s->buf;
        int64_t k = s->patch_children * (cbefore + 1);
        for (c = 0; c < s->cells_per_patch; k += s->patch_children, ++c)
        {
            *idata++ = k;
        }
    }
    write_buffer (s, s->psize_offsets);
}

static void
write_types_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int c;

    char *cdata = s->buf;
    if(s->dim == 2)
    {
        for (c = 0; c < s->cells_per_patch; ++c)
        {
            *cdata++ = 9;
        }
    }
    else 
    {
        for (c = 0; c < s->cells_per_patch; ++c)
        {
            *cdata++ = 12;
        }
    }
    write_buffer (s, s->psize_types);
}

static void
write_mpirank_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                  int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int c;

    int *idata = (int *) s->buf;
    for (c = 0; c < s->cells_per_patch; ++c)
    {
        *idata++ = domain->mpirank;
    }
    write_buffer (s, s->psize_mpirank);
}

static void
write_blockno_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                  int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int c;

    int *idata = (int *) s->buf;
    for (c = 0; c < s->cells_per_patch; ++c)
    {
        *idata++ = blockno;
    }
    write_buffer (s, s->psize_blockno);
}

static void
write_patchno_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                  int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;
    int c;
    const int64_t gpno =
        domain->global_num_patches_before +
        domain->blocks[blockno].num_patches_before + patchno;

    if (s->fits32)
    {
        const int32_t igpno = (int32_t) gpno;
        int32_t *idata = (int32_t *) s->buf;
        for (c = 0; c < s->cells_per_patch; ++c)
        {
            *idata++ = igpno;
        }
    }
    else
    {
        int64_t *idata = (int64_t *) s->buf;
        for (c = 0; c < s->cells_per_patch; ++c)
        {
            *idata++ = gpno;
        }
    }
    write_buffer (s, s->psize_patchno);
}

static void
write_3d_patch_points (fclaw_global_t * glob,
                       fclaw_patch_t * patch,
                       int blockno,
                       int patchno,
                       double * points)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fclaw2d_map_context_t* cont = fclaw2d_map_get(glob);
    /* Enumerate point coordinates in the patch */
    int i, j, k;
    double xpp,ypp,zpp;
    for (k = 0; k <= mz; ++k)
    {
        const double z = zlower + k * dz;
        for (j = 0; j <= my; ++j)
        {
            const double y = ylower + j * dy;
            for (i = 0; i <= mx; ++i)
            {
                const double x = xlower + i * dx;
                if (fclaw_opt->manifold)
                {
                    FCLAW3D_MAP_C2M(&cont,&blockno,&x,&y,&z,&xpp,&ypp,&zpp);
                    *points++ = xpp;
                    *points++ = ypp;
                    *points++ = zpp;
                }
                else
                {
                    *points++ = x;
                    *points++ = y;
                    *points++ = z;
                }
            }
        }
    }
}

static void
write_patch_points (fclaw_global_t * glob,
                    fclaw_patch_t * patch,
                    int blockno,
                    int patchno,
                    double * points)
{
    write_3d_patch_points(glob, patch, blockno, patchno, points);
}

static void
write_3d_patch_connectivity (fclaw_global_t * glob,
                             fclaw_patch_t * patch,
                             int blockno,
                             int patchno,
                             long num_points_before,
                             long * conn)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);
    for (int k = 0; k < mz; ++k)
    {
        for (int j = 0; j < my; ++j)
        {
            for (int i = 0; i < mx; ++i)
            {
                long l = num_points_before + i + j * (mx + 1)
                     + k * (my + 1) * (mx + 1);
                *conn++ = l;
                *conn++ = l + 1;
                *conn++ = l + (mx + 2);
                *conn++ = l + (mx + 1);
                *conn++ = l + (mx + 1) * (my + 1);
                *conn++ = l + (mx + 1) * (my + 1) + 1;
                *conn++ = l + (mx + 1) * (my + 1) + (mx + 2);
                *conn++ = l + (mx + 1) * (my + 1) + (mx + 1);
            }
        }
    }
}
static void
write_patch_q (fclaw_global_t * glob,
               fclaw_patch_t * patch,
               int blockno,
               int patchno,
               double * q_out)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    int meqn;
    double* q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    for (int k = 0; k < mz; ++k)
    {
        for (int j = 0; j < my; ++j)
        {
            for (int i = 0; i < mx; ++i)
            {
                for (int m = 0; m < meqn; m++)
                {
                    size_t out_idx = k*mx*my*meqn + j*mx*meqn + i*meqn + m;
                    size_t in_idx = m*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc) + (k+mbc)*(mx+2*mbc)*(my+2*mbc) + (j+mbc)*(mx+2*mbc) + (i+mbc);
                    q_out[out_idx] = q[in_idx];
                }
            }
        }
    }
}
static void
write_patch_connectivity (fclaw_global_t * glob,
                          fclaw_patch_t * patch,
                          int blockno,
                          int patchno,
                          long num_points_before,
                          long * connectivity)
{
    write_3d_patch_connectivity(glob, patch, blockno, patchno, num_points_before, connectivity);
}

static
herr_t
set_attribute_numerical(hid_t loc_id, const char *obj_name, const char *attr_name, size_t size,
                        hid_t tid, const void *data)
{

    hid_t   obj_id, sid, attr_id;
    hsize_t dim_size = size;
    htri_t  has_attr;

    /* check the arguments */
    if (obj_name == NULL)
        return -1;
    if (attr_name == NULL)
        return -1;

    /* Open the object */
    if ((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
        return -1;

    /* Create the data space for the attribute. */
    if ((sid = H5Screate_simple(1, &dim_size, NULL)) < 0)
        goto out;

    /* Delete the attribute if it already exists */
    if ((has_attr = H5Aexists(obj_id, attr_name)) < 0)
        goto out;
    if (has_attr > 0)
        if (H5Adelete(obj_id, attr_name) < 0)
            goto out;

    /* Create the attribute. */
    if ((attr_id = H5Acreate2(obj_id, attr_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        goto out;

    /* Write the attribute data. */
    if (H5Awrite(attr_id, tid, data) < 0)
        goto out;

    /* Close the attribute. */
    if (H5Aclose(attr_id) < 0)
        goto out;

    /* Close the dataspace. */
    if (H5Sclose(sid) < 0)
        goto out;

    /* Close the object */
    if (H5Oclose(obj_id) < 0)
        return -1;

    return 0;

out:
    H5Oclose(obj_id);
    return -1;
}
static herr_t
make_single_value_dataset_numerical(int mpirank,
                                    hid_t loc_id, 
                                    const char *dset_name, 
                                    int rank, 
                                    const hsize_t *dims, 
                                    hid_t tid,
                                    const void *data)
{
    hid_t did = -1, sid = -1;

    /* check the arguments */
    if (dset_name == NULL)
        return -1;

    /* Create the data space for the dataset. */
    if ((sid = H5Screate_simple(rank, dims, NULL)) < 0)
        return -1;

    /* Create the dataset. */
    if ((did = H5Dcreate2(loc_id, dset_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        goto out;

    /* Write the dataset only if there is data to write */
    if (data && mpirank == 0)
        if (H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) < 0)
            goto out;

    /* End access to the dataset and release resources used by it. */
    if (H5Dclose(did) < 0)
        return -1;

    /* Terminate access to the data space. */
    if (H5Sclose(sid) < 0)
        return -1;

    return 0;

out:
    H5E_BEGIN_TRY
    {
        H5Dclose(did);
        H5Sclose(sid);
    }
    H5E_END_TRY
    return -1;
}
static herr_t
make_dataset_numerical(hid_t loc_id, 
                       const char *dset_name, 
                       int rank, 
                       const hsize_t *dims, 
                       const hsize_t *chunk_dims,
                       const hsize_t* slab_start, 
                       const hsize_t* slab_dims, 
                       hid_t tid,
                       const void *data)
{
    hid_t did = -1, sid = -1, prop_id = -1, plist_id = -1;

    /* check the arguments */
    if (dset_name == NULL)
        return -1;
    if((prop_id = H5Pcreate(H5P_DATASET_CREATE)) < 0)
        return -1;
    
    if(chunk_dims != NULL)
    {
        if(H5Pset_chunk(prop_id, rank, chunk_dims) < 0)
            return -1;
        if(H5Pset_deflate(prop_id, 5) < 0)
            return -1;
    }

    
    /* Create the data space for the dataset. */
    if ((sid = H5Screate_simple(rank, dims, NULL)) < 0)
        return -1;


    /* Create the dataset. */
    if ((did = H5Dcreate2(loc_id, dset_name, tid, sid, H5P_DEFAULT, prop_id, H5P_DEFAULT)) < 0)
        goto out;

    hid_t memspace = H5Screate_simple(rank, slab_dims, NULL);
    if(memspace < 0)
        return -1;

    hid_t filespace = H5Dget_space(did);
    if(filespace < 0)
        return -1;

    if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET, slab_start, NULL, slab_dims, NULL) < 0)
        return -1;

    if((plist_id = H5Pcreate(H5P_DATASET_XFER)) < 0)
        return -1;

    if(H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE) < 0)
       return -1;

    /* Write the dataset only if there is data to write */
    if (data)
        if (H5Dwrite(did, tid, memspace, filespace, plist_id, data) < 0)
            goto out;

    /* End access to the dataset and release resources used by it. */
    if (H5Dclose(did) < 0)
        return -1;

    /* Terminate access to the data space. */
    if (H5Sclose(sid) < 0)
        return -1;

    if (H5Pclose(prop_id) < 0)
        return -1;

    return 0;

out:
    H5E_BEGIN_TRY
    {
        H5Dclose(did);
        H5Sclose(sid);
        H5Sclose(prop_id);
    }
    H5E_END_TRY
    return -1;
}
//this is copied form hdf5lt needed to change to NULLPAD for VTK
herr_t
set_attribute_string(hid_t loc_id, const char *obj_name, const char *attr_name, const char *attr_data)
{
    hid_t  attr_type;
    hid_t  attr_space_id;
    hid_t  attr_id;
    hid_t  obj_id;
    htri_t has_attr;
    size_t attr_size;

    /* check the arguments */
    if (obj_name == NULL)
        return -1;
    if (attr_name == NULL)
        return -1;
    if (attr_data == NULL)
        return -1;

    /* Open the object */
    if ((obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT)) < 0)
        return -1;

    /* Create the attribute */
    if ((attr_type = H5Tcopy(H5T_C_S1)) < 0)
        goto out;

    attr_size = strlen(attr_data);

    if (H5Tset_size(attr_type, (size_t)attr_size) < 0)
        goto out;

    if (H5Tset_strpad(attr_type, H5T_STR_NULLPAD) < 0)
        goto out;

    if ((attr_space_id = H5Screate(H5S_SCALAR)) < 0)
        goto out;

    /* Delete the attribute if it already exists */
    if ((has_attr = H5Aexists(obj_id, attr_name)) < 0)
        goto out;
    if (has_attr > 0)
        if (H5Adelete(obj_id, attr_name) < 0)
            goto out;

    /* Create and write the attribute */

    if ((attr_id = H5Acreate2(obj_id, attr_name, attr_type, attr_space_id, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        goto out;

    if (H5Awrite(attr_id, attr_type, attr_data) < 0)
        goto out;

    if (H5Aclose(attr_id) < 0)
        goto out;

    if (H5Sclose(attr_space_id) < 0)
        goto out;

    if (H5Tclose(attr_type) < 0)
        goto out;

    /* Close the object */
    if (H5Oclose(obj_id) < 0)
        return -1;

    return 0;

out:

    H5Oclose(obj_id);
    return -1;
}
static int
fclaw_hdf5_write_file (int dim, fclaw_global_t * glob, const char *basename,
                      int mx, int my, int mz,
                      int meqn)
{
    const fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int ierr;

    char vtkhdf[8] = "/VTKHDF";
    char celldata[18] = "/VTKHDF/CellData";
    
    // Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if(H5Pset_fapl_mpio(plist_id, glob->mpicomm, MPI_INFO_NULL) < 0)
        return -1;

    // Create a new file collectively and release property list identifier.
    hid_t file_id = H5Fcreate(basename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    if(file_id < 0)
        return -1;

    H5Pclose(plist_id);

    hid_t gid1 = H5Gcreate2(file_id, vtkhdf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    int num_cells_per_patch = mx * my * mz;
    int num_points_per_patch = (mx + 1) * (my + 1) * (mz + 1);
    int local_num_patches = glob->domain->local_num_patches;
    int global_num_patches = glob->domain->global_num_patches;

    
    int num_points_per_cell;
    uint8_t cell_type;
    if(clawpatch_opt->patch_dim == 2)
    {
        num_points_per_cell = 4;
        cell_type = 9;
    }
    else
    {
        num_points_per_cell = 8;
        cell_type = 12; 
    }

    int vtk_version[2] = {1, 0};
    set_attribute_numerical(file_id, vtkhdf, "Version", 2, H5T_NATIVE_INT, vtk_version);
    //hid_t attr_id = H5Acreate (gid1, "Units", H5T_STD_I32BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    //H5Aclose (attr_id);
    set_attribute_string(file_id, vtkhdf, "Type", "UnstructuredGrid");
    //const char* type = "UnstructuredGrid";
    //H5LTset_attribute(fid, vtkhdf, "Type", type);
    
    // The idea here is to write a paritioned vtu, where each partition is a patch
    // this will make processing easier in matlab

    hsize_t dims[3] = {0,0,0};

    // write single value datasets for vtk

    long number_of_cells = global_num_patches * num_cells_per_patch;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, gid1, "NumberOfCells", 1, dims, H5T_NATIVE_LONG, &number_of_cells);

    long number_of_points = global_num_patches * num_points_per_patch;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, gid1, "NumberOfPoints", 1, dims, H5T_NATIVE_LONG, &number_of_points);


    long number_of_connectivity_ids = global_num_patches * num_cells_per_patch * num_points_per_cell;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, gid1, "NumberOfConnectivityIds", 1, dims, H5T_NATIVE_LONG, &number_of_connectivity_ids);



    hsize_t chunk_dims[3] = {0,0,0};
    hsize_t slab_dims[3] = {0,0,0};
    hsize_t slab_start[3] = {0,0,0};

    uint8_t types[local_num_patches * num_cells_per_patch]; 
    for(int i = 0; i < local_num_patches * num_cells_per_patch; i++){
        types[i] = cell_type;
    }

    dims[0] = global_num_patches * num_cells_per_patch;
    chunk_dims[0] = num_cells_per_patch;
    slab_start[0] = glob->domain->global_num_patches_before * num_cells_per_patch;
    slab_dims[0] = local_num_patches * num_cells_per_patch;
    make_dataset_numerical(gid1, "Types", 1, dims, chunk_dims, slab_start, slab_dims, H5T_NATIVE_UINT8, types);


#if 0

    // write offsets
    // calculate additional end lenth for rank == size -1
    int end = (glob->mpirank == glob->mpisize - 1) ? 1 : 0;
    long *offsets = FCLAW_ALLOC(long,local_num_patches * num_cells_per_patch + end);
    int curr_offset = glob->domain->global_num_patches_before * num_cells_per_patch * num_points_per_cell;
    for(int i=0; i < local_num_patches * num_cells_per_patch + end; i++){
        offsets[i] = curr_offset;
        curr_offset += num_points_per_cell;
    }

    dims[0] = global_num_patches * num_cells_per_patch + 1;
    chunk_dims[0] = num_cells_per_patch;
    slab_start[0] = glob->domain->global_num_patches_before * num_cells_per_patch;
    slab_dims[0] = local_num_patches * num_cells_per_patch + end;
    make_dataset_numerical(gid1, "Offsets", 1, dims, chunk_dims, slab_start, slab_dims, H5T_NATIVE_LONG, offsets);
    FCLAW_FREE(offsets);


    double* points = FCLAW_ALLOC(double, global_num_patches * num_points_per_patch * 3);
    dims[0] = global_num_patches * num_points_per_patch;
    dims[1] = 3;
    chunk_dims[0] = num_points_per_patch;
    chunk_dims[1] = 3;
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++)
    {
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++)
        {
            fclaw_patch_t* patch = &block->patches[patchno];
            double *patch_points = &points[(block->num_patches_before + patchno) * num_points_per_patch * 3];
            write_patch_points(glob, patch, blockno, patchno, patch_points);
        }
    }

    make_dataset_numerical(gid1, "Points", 2, dims, chunk_dims, H5T_NATIVE_DOUBLE, points);
    FCLAW_FREE(points);

    long* connectivity = FCLAW_ALLOC(long, global_num_patches * num_cells_per_patch * num_points_per_cell);

    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++)
    {
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++)
        {
            fclaw_patch_t* patch = &block->patches[patchno];
            long num_points_before = (block->num_patches_before + patchno) * num_points_per_patch;
            long *patch_connectivity = &connectivity[(block->num_patches_before + patchno) * num_cells_per_patch * num_points_per_cell];
            write_patch_connectivity(glob, patch, blockno, patchno, num_points_before, patch_connectivity);
        }
    }

    dims[0] = global_num_patches * num_cells_per_patch * num_points_per_cell;
    chunk_dims[0] = num_cells_per_patch * num_points_per_cell;
    make_dataset_numerical(gid1, "Connectivity", 1, dims, chunk_dims, H5T_NATIVE_LONG, connectivity);
    FCLAW_FREE(connectivity);

    /* avoid resource leaks by closing */
    H5Gclose(gid1);




    gid1 = H5Gcreate2(file_id, celldata, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    //create contiguous array for q
    double * q = FCLAW_ALLOC(double, global_num_patches * num_cells_per_patch * meqn);
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++)
    {
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++)
        {
            fclaw_patch_t* patch = &block->patches[patchno];
            double *patch_q = &q[(block->num_patches_before + patchno) * num_cells_per_patch * meqn];
            write_patch_q(glob, patch, blockno, patchno, patch_q);
        }
    }
    dims[0] = num_cells_per_patch*global_num_patches;
    dims[1] = meqn;
    chunk_dims[0] = num_cells_per_patch;
    chunk_dims[1] = meqn;
    make_dataset_numerical(gid1, "meqn", 1, dims, chunk_dims, H5T_NATIVE_DOUBLE, q);
    FCLAW_FREE(q);

    // write blockno
    int *blockno_array = FCLAW_ALLOC(int, global_num_patches * num_cells_per_patch);
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++)
    {
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++)
        {
            int *patch_blockno = &blockno_array[(block->num_patches_before + patchno) * num_cells_per_patch];
            for(int i = 0; i < num_cells_per_patch; i++){
                patch_blockno[i] = blockno;
            }
        }
    }
    dims[0] = global_num_patches * num_cells_per_patch;
    chunk_dims[0] = num_cells_per_patch;
    make_dataset_numerical(gid1, "blockno", 1, dims, chunk_dims, H5T_NATIVE_INT, blockno_array);
    FCLAW_FREE(blockno_array);

    //write patchno
    int *patchno_array = FCLAW_ALLOC(int, global_num_patches * num_cells_per_patch);
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++)
    {
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++)
        {
            int *patch_patchno = &patchno_array[(block->num_patches_before + patchno) * num_cells_per_patch];
            for(int i = 0; i < num_cells_per_patch; i++){
                patch_patchno[i] = patchno;
            }
        }
    }
    dims[0] = global_num_patches * num_cells_per_patch;
    chunk_dims[0] = num_cells_per_patch;
    make_dataset_numerical(gid1, "patchno", 1, dims, chunk_dims, H5T_NATIVE_INT, patchno_array);

    //write mpirank
    int *mpirank_array = FCLAW_ALLOC(int, global_num_patches * num_cells_per_patch);
    for(int blockno = 0; blockno < glob->domain->num_blocks; blockno++)
    {
        fclaw_block_t* block = &glob->domain->blocks[blockno];
        for(int patchno = 0; patchno < block->num_patches; patchno++)
        {
            int *patch_mpirank = &mpirank_array[(block->num_patches_before + patchno) * num_cells_per_patch];
            for(int i = 0; i < num_cells_per_patch; i++){
                patch_mpirank[i] = glob->mpirank;
            }
        }
    }
    dims[0] = global_num_patches * num_cells_per_patch;
    chunk_dims[0] = num_cells_per_patch;
    make_dataset_numerical(gid1, "mpirank", 1, dims, chunk_dims, H5T_NATIVE_INT, mpirank_array);
    
#endif
    H5Gclose(gid1);

    H5Fclose(file_id);
    
    return EXIT_SUCCESS;
}

int
fclaw_hdf5_write_2d_file (fclaw_global_t * glob, const char *basename,
                        int mx, int my,
                        int meqn)
{
    return fclaw_hdf5_write_file(2,glob,basename,mx,my,0,meqn);
}

int
fclaw_hdf5_write_3d_file (fclaw_global_t * glob, const char *basename,
                        int mx, int my, int mz,
                        int meqn)
{
    return fclaw_hdf5_write_file(3,glob,basename,mx,my,mz,meqn);
}

/*  ---------------------------------------------------------------------------
    Public interface
    --------------------------------------------------------------------------- */

void fclaw_clawpatch_output_hdf5_to_file (fclaw_global_t * glob, const char* filename)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);


    if(clawpatch_opt->patch_dim == 2)
    {
        fclaw_hdf5_write_2d_file (glob, filename,
                                clawpatch_opt->mx, 
                                clawpatch_opt->my,
                                clawpatch_opt->meqn);

    }
    else 
    {
        fclaw_hdf5_write_3d_file (glob, filename,
                                clawpatch_opt->mx, 
                                clawpatch_opt->my, 
                                clawpatch_opt->mz,
                                clawpatch_opt->meqn);
    }
}
void fclaw_clawpatch_output_hdf5 (fclaw_global_t * glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d.hdf", fclaw_opt->prefix, iframe);

    fclaw_clawpatch_output_hdf5_to_file(glob,basename);
}



