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

static int
fclaw2d_hdf5_write_header (fclaw_domain_t * domain, fclaw2d_hdf5_state_t * s)
{
    int retval;
    FILE *file;

    /* unconditionally open the file */
    file = fopen (s->filename, "wb");
    if (file == NULL)
    {
        return -1;
    }
    s->file = file;

    /* stop writing after first unsuccessful operation */
    retval = 0;
    retval = retval || fprintf (file, "<?xml version=\"1.0\"?>\n") < 0;
    retval = retval || fprintf (file, "<VTKFile type=\"UnstructuredGrid\" "
                                "version=\"0.1\" "
                                "byte_order=\"LittleEndian\" "
                                "header_type=\"UInt64\" "
                                ">\n") < 0;
    retval = retval || fprintf (file, " <UnstructuredGrid>\n") < 0;
    retval = retval || fprintf (file, "  <Piece NumberOfPoints=\"%lld\" "
                                "NumberOfCells=\"%lld\">\n",
                                (long long) s->global_num_points,
                                (long long) s->global_num_cells) < 0;
    retval = retval || fprintf (file, "   <Points>\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"Float64\" "
                                "Name=\"Position\" NumberOfComponents=\"3\" "
                                "format=\"appended\" offset=\"%lld\">\n",
                                (long long) s->offset_position) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "   </Points>\n") < 0;
    retval = retval || fprintf (file, "   <Cells>\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"%s\" "
                                "Name=\"connectivity\" format=\"appended\" "
                                "offset=\"%lld\">\n", s->inttype,
                                (long long) s->offset_connectivity) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"%s\" "
                                "Name=\"offsets\" format=\"appended\" "
                                "offset=\"%lld\">\n", s->inttype,
                                (long long) s->offset_offsets) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"UInt8\" "
                                "Name=\"types\" format=\"appended\" "
                                "offset=\"%lld\">\n",
                                (long long) s->offset_types) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "   </Cells>\n") < 0;
    retval = retval || fprintf (file, "   <CellData Scalars=\"mpirank,"
                                "blockno,patchno\" Fields=\"meqn\">\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"Int32\" "
                                "Name=\"mpirank\" format=\"appended\" "
                                "offset=\"%lld\">\n",
                                (long long) s->offset_mpirank) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"Int32\" "
                                "Name=\"blockno\" format=\"appended\" "
                                "offset=\"%lld\">\n",
                                (long long) s->offset_blockno) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "    <DataArray type=\"%s\" "
                                "Name=\"patchno\" format=\"appended\" "
                                "offset=\"%lld\">\n", s->inttype,
                                (long long) s->offset_patchno) < 0;
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    if (s->meqn == 1)
    {
        retval = retval || fprintf (file, "    <DataArray type=\"Float32\" "
                                    "Name=\"meqn\" "
                                    "format=\"appended\" offset=\"%lld\">\n",
                                    (long long) s->offset_meqn) < 0;
    }
    else
    {
        retval = retval || fprintf (file, "    <DataArray type=\"Float32\" "
                                    "Name=\"meqn\" NumberOfComponents=\"%d\" "
                                    "format=\"appended\" offset=\"%lld\">\n",
                                    s->meqn, (long long) s->offset_meqn) < 0;
    }
    retval = retval || fprintf (file, "    </DataArray>\n") < 0;
    retval = retval || fprintf (file, "   </CellData>\n") < 0;
    retval = retval || fprintf (file, "   <PointData>\n") < 0;
    retval = retval || fprintf (file, "   </PointData>\n") < 0;
    retval = retval || fprintf (file, "  </Piece>\n") < 0;
    retval = retval || fprintf (file, " </UnstructuredGrid>\n") < 0;
    retval = retval || fprintf (file, " <AppendedData "
                                "encoding=\"raw\">\n  _") < 0;

#ifdef P4EST_ENABLE_MPIIO
    /* unconditionally close the file when in MPI I/O mode */
    retval = fclose (file) || retval;
    s->file = NULL;
#endif

    return retval ? -1 : 0;
}

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
write_meqn_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
               int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;

    fclaw2d_hdf5_state_t *s = (fclaw2d_hdf5_state_t *) g->user;

    s->value_cb (g->glob, patch, blockno, patchno, s->buf);
    write_buffer (s, s->psize_meqn);
}

static void
fclaw2d_hdf5_write_field (fclaw_global_t * glob, fclaw2d_hdf5_state_t * s,
                         int64_t offset_field, int64_t psize_field,
                         fclaw_patch_callback_t cb)
{
    fclaw_domain_t *domain = glob->domain;

    int64_t bcount;
#ifndef P4EST_ENABLE_MPIIO
    size_t retvalz;
#else
    int mpiret;
    MPI_Offset mpipos;
#ifdef P4EST_ENABLE_DEBUG
    MPI_Offset mpinew;
#endif
    MPI_Status mpistatus;
#endif

    s->buf = P4EST_ALLOC (char, psize_field);
#ifdef P4EST_ENABLE_MPIIO
    mpipos = s->mpibegin + offset_field;
    if (domain->mpirank > 0)
    {
        /* account for byte count */
        mpipos += s->ndsize + psize_field * domain->global_num_patches_before;
    }
    mpiret = MPI_File_seek (s->mpifile, mpipos, MPI_SEEK_SET);
    SC_CHECK_MPI (mpiret);
#endif
    if (domain->mpirank == 0)
    {
        /* write byte count */
        bcount = psize_field * domain->global_num_patches;
#if 0
        P4EST_LDEBUGF ("offset %lld psize %lld bcount %d %o %x\n",
                       (long long) offset_field, (long long) psize_field,
                       bcount, bcount, bcount);
#endif
#ifndef P4EST_ENABLE_MPIIO
        retvalz = fwrite (&bcount, s->ndsize, 1, s->file);
        SC_CHECK_ABORT (retvalz == 1, "VTK file write failed");
#else
        mpiret = MPI_File_write (s->mpifile, &bcount, 1, MPI_LONG, &mpistatus);
        SC_CHECK_MPI (mpiret);
#endif
    }
    fclaw_global_iterate_patches (glob, cb, s);
    P4EST_FREE (s->buf);

#ifdef P4EST_ENABLE_MPIIO
#ifdef P4EST_ENABLE_DEBUG
    mpiret = MPI_File_get_position (s->mpifile, &mpinew);
    SC_CHECK_MPI (mpiret);
    P4EST_ASSERT (mpinew - mpipos ==
                  (domain->mpirank == 0 ? s->ndsize : 0) +
                  domain->local_num_patches * psize_field);
    P4EST_ASSERT (domain->mpirank < domain->mpisize - 1 ||
                  mpinew - s->mpibegin ==
                  offset_field + s->ndsize +
                  psize_field * domain->global_num_patches);
#endif
#endif
}

static void
fclaw2d_hdf5_write_data (fclaw_global_t * glob, fclaw2d_hdf5_state_t * s)
{
#ifdef P4EST_ENABLE_MPIIO
    int mpiret;
    MPI_Offset mpipos;

    /* collectively open the file in append mode and reserve space */
    mpiret = MPI_File_open (glob->mpicomm, s->filename,
                            MPI_MODE_RDWR | MPI_MODE_APPEND |
                            MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &s->mpifile);
    SC_CHECK_MPI (mpiret);
    mpiret = MPI_File_get_position (s->mpifile, &s->mpibegin);
    SC_CHECK_MPI (mpiret);
    mpipos = s->mpibegin + (MPI_Offset) s->offset_end;
    mpiret = MPI_File_preallocate (s->mpifile, mpipos);
    SC_CHECK_MPI (mpiret);
#endif

    /* write meta data fields */
    fclaw2d_hdf5_write_field (glob, s, s->offset_position, s->psize_position,
                             write_position_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_connectivity,
                             s->psize_connectivity, (s->dim == 2) ? write_2d_connectivity_cb : write_3d_connectivity_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_offsets, s->psize_offsets,
                             write_offsets_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_types, s->psize_types,
                             write_types_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_mpirank, s->psize_mpirank,
                             write_mpirank_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_blockno, s->psize_blockno,
                             write_blockno_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_patchno, s->psize_patchno,
                             write_patchno_cb);
    fclaw2d_hdf5_write_field (glob, s, s->offset_meqn, s->psize_meqn,
                             write_meqn_cb);

#ifdef P4EST_ENABLE_MPIIO
    /* collectively close the file */
    mpiret = MPI_File_close (&s->mpifile);
    SC_CHECK_MPI (mpiret);
#endif
}

static int
fclaw2d_hdf5_write_footer (fclaw_domain_t * domain, fclaw2d_hdf5_state_t * s)
{
    int retval;
    FILE *file;

#ifndef P4EST_ENABLE_MPIIO
    file = s->file;
#else
    /* unconditionally open the file */
    file = fopen (s->filename, "ab");
    if (file == NULL)
    {
        return -1;
    }
    s->file = file;
#endif

    /* stop writing after first unsuccessful operation */
    retval = 0;
    retval = retval || fprintf (file, "\n </AppendedData>\n</VTKFile>\n") < 0;

    /* unconditionally close the file */
    retval = fclose (file) || retval;
    s->file = NULL;

    return retval ? -1 : 0;
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

static int
fclaw_hdf5_write_file (int dim, fclaw_global_t * glob, const char *basename,
                      int mx, int my, int mz,
                      int meqn,
                      double hdf5space, int hdf5write,
                      fclaw_hdf5_patch_data_t coordinate_cb,
                      fclaw_hdf5_patch_data_t value_cb)
{
    const fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int ierr;

    char vtkhdf[8] = "/VTKHDF";
    char pointdata[18] = "/VTKHDF/PointData";
    
    /* create HDF5 file */
    hid_t fid = H5Fcreate(basename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    hid_t gid1 = H5Gcreate2(fid, vtkhdf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
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
    H5LTset_attribute_int(fid, vtkhdf, "Version", vtk_version, 2);
    
    // The idea here is to write a paritioned vtu, where each partition is a patch
    // this will make processing easier in matlab

    int number_of_cells[global_num_patches];
    for(int i = 0; i < global_num_patches; i++){
        number_of_cells[i] = num_cells_per_patch;
    }

    hsize_t dims[3] = {0,0,0};
    dims[0] = global_num_patches;
    H5LTmake_dataset_int(gid1, "NumberOfCells", 1, dims, number_of_cells);

    int number_of_points[global_num_patches];
    for(int i = 0; i < global_num_patches; i++){
        number_of_points[i] = num_points_per_patch;
    }

    dims[0] = global_num_patches;
    H5LTmake_dataset_int(gid1, "NumberOfPoints", 1, dims, number_of_points);


    int number_of_connectivity_ids[global_num_patches];
    for(int i = 0; i < global_num_patches; i++){
        number_of_connectivity_ids[i] = num_cells_per_patch * num_points_per_cell;
    }

    dims[0] = global_num_patches;
    H5LTmake_dataset_int(gid1, "NumberOfConnectivityIds", 1, dims, number_of_connectivity_ids);


    char types[global_num_patches * num_cells_per_patch]; 
    for(int i = 0; i < global_num_patches * num_cells_per_patch; i++){
        types[i] = cell_type;
    }

    dims[0] = global_num_patches * num_cells_per_patch;
    H5LTmake_dataset_char(gid1, "Types", 1, dims, types);


    int offsets[global_num_patches * num_cells_per_patch + 1];
    int curr_offset = 0;
    for(int i = 0; i < global_num_patches * num_cells_per_patch + 1; i++){
        offsets[i] = curr_offset;
        curr_offset += num_points_per_cell;
    }

    dims[0] = global_num_patches * num_cells_per_patch + 1;
    H5LTmake_dataset_int(gid1, "Offsets", 1, dims, offsets);


    double* points = FCLAW_ALLOC(double, global_num_patches * num_points_per_patch * 3);
    hsize_t points_dim[3] = {global_num_patches, num_points_per_patch, 3};
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

    H5LTmake_dataset_double(gid1, "Points", 3, points_dim, points);
    FCLAW_FREE(points);

    /* avoid resource leaks by closing */
    H5Gclose(gid1);




    //H5LTset_attribute_int(fid, vtkhdf, "NumberOfConnectivityIds", number_of_points, global_num_patches);

    //float spacing[3] = {0.131579, 0.125, 0.0952381};
    //H5LTset_attribute_float(fid, vtkhdf, "Spacing", spacing, 3);
    
    //int whole_extent[6] = {0, 19, 0, 20, 0, 21};
    //H5LTset_attribute_int(fid, vtkhdf, "WholeExtent", whole_extent, 6);
    
    //H5LTset_attribute_string(fid, vtkhdf, "Type", "UnstructuredGrid");
    
    
    //gid1 = H5Gcreate2(fid, pointdata, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    //H5LTset_attribute_string(fid, pointdata, "Scalars", "Iterations");
    
    ///* declare 3D array of test data */
    //hsize_t dims[3] = {20, 21, 22};
    ///* malloc array */
    //float* Iterations = (float*) malloc(dims[0] * dims[1] * dims[2] * sizeof(float));
    
    ///* arbitrary array values */
    //for(int i = 0; i < dims[0] * dims[1] * dims[2]; i++){
    //  Iterations[i] = (float) rand() / (float) RAND_MAX;
    //}
    
    //H5LTmake_dataset_float(gid1, "Iterations", 3, dims, Iterations);
    
    //H5Gclose(gid1);
    

    H5Fclose(fid);
    
    return EXIT_SUCCESS;
}

int
fclaw_hdf5_write_2d_file (fclaw_global_t * glob, const char *basename,
                        int mx, int my,
                        int meqn,
                        double hdf5space, int hdf5write,
                        fclaw_hdf5_patch_data_t coordinate_cb,
                        fclaw_hdf5_patch_data_t value_cb)
{
    return fclaw_hdf5_write_file(2,glob,basename,mx,my,0,meqn,hdf5space,hdf5write,
                                coordinate_cb,value_cb);
}

int
fclaw_hdf5_write_3d_file (fclaw_global_t * glob, const char *basename,
                        int mx, int my, int mz,
                        int meqn,
                        double hdf5space, int hdf5write,
                        fclaw_hdf5_patch_data_t coordinate_cb,
                        fclaw_hdf5_patch_data_t value_cb)
{
    return fclaw_hdf5_write_file(3,glob,basename,mx,my,mz,meqn,hdf5space,hdf5write,
                                coordinate_cb,value_cb);
}

static void
fclaw2d_output_hdf5_coordinate_cb (fclaw_global_t * glob,
                                  fclaw_patch_t * patch,
                                  int blockno, int patchno,
                                  char *a)
{
    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fclaw2d_map_context_t* cont = fclaw2d_map_get(glob);

    /* Enumerate point coordinates in the patch */
    double *d = (double *) a;
    int i, j;
    double xpp,ypp,zpp;
    for (j = 0; j <= my; ++j)
    {
        const double y = ylower + j * dy;
        for (i = 0; i <= mx; ++i)
        {
            const double x = xlower + i * dx;
            if (fclaw_opt->manifold)
            {
                FCLAW2D_MAP_C2M(&cont,&blockno,&x,&y,&xpp,&ypp,&zpp);
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
fclaw3d_output_hdf5_coordinate_cb (fclaw_global_t * glob,
                                  fclaw_patch_t * patch,
                                  int blockno, int patchno,
                                  char *a)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fclaw2d_map_context_t* cont = fclaw2d_map_get(glob);
    /* Enumerate point coordinates in the patch */
    double *d = (double *) a;
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
fclaw2d_output_hdf5_value_cb (fclaw_global_t * glob,
                             fclaw_patch_t * patch,
                             int blockno, int patchno,
                             char *a)
{
    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;

    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const int xlane = mx + 2 * mbc;
    const int ylane = my + 2 * mbc;

    // Enumerate equation data in the patch
    float *f = (float *) a;
    int i, j, k;
    for (j = 0; j < my; ++j)
    {
        for (i = 0; i < mx; ++i)
        {
            for (k = 0; k < meqn; ++k)
            {
                /* For Clawpack 5.0 layout */
                //*f++ = (float) q[((j+mbc)*xlane + (i+mbc))*meqn + k];

                /* For Clawpack 4.x layout */
                *f++ = (float) q[(k * ylane + j + mbc) * xlane + i + mbc];
            }
        }
    }
}

static void
fclaw3d_output_hdf5_value_cb (fclaw_global_t * glob,
                             fclaw_patch_t * patch,
                             int blockno, int patchno,
                             char *a)
{

    int meqn;
    double *q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    int mx,my,mz,mbc;
    double xlower,ylower,zlower,dx,dy,dz;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                               &xlower,&ylower,&zlower, &dx,&dy, &dz);

    const int xlane = mx + 2 * mbc;
    const int ylane = my + 2 * mbc;
    const int zlane = mz + 2 * mbc;

    // Enumerate equation data in the patch
    float *f = (float *) a;
    int i, j, k, eqn;
    for (k = 0; k < mz; ++k)
    {
        for (j = 0; j < my; ++j)
        {
            for (i = 0; i < mx; ++i)
            {
                for (eqn = 0; eqn < meqn; ++eqn)
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

/*  --------------------------------------------------------------------------
    Used for debugging
    ------------------------------------------------------------------------- */
#if 0
static void
fclaw2d_output_write_hdf5_debug (fclaw_global_t * glob, const char *basename)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    (void) fclaw2d_hdf5_write_file (glob, basename,
                                   clawpatch_opt->mx, clawpatch_opt->my,
#if PATCH_DIM == 3
                                   clawpatch_opt->mz,
#endif
                                   clawpatch_opt->meqn,
                                   fclaw_opt->hdf5space, 0,
                                   fclaw2d_output_hdf5_coordinate_cb,
                                   fclaw2d_output_hdf5_value_cb);
}
#endif


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
                                clawpatch_opt->meqn,
                                fclaw_opt->vtkspace, 0,
                                fclaw2d_output_hdf5_coordinate_cb,
                                fclaw2d_output_hdf5_value_cb);

    }
    else 
    {
        fclaw_hdf5_write_3d_file (glob, filename,
                                clawpatch_opt->mx, 
                                clawpatch_opt->my, 
                                clawpatch_opt->mz,
                                clawpatch_opt->meqn,
                                fclaw_opt->vtkspace, 0,
                                fclaw3d_output_hdf5_coordinate_cb,
                                fclaw3d_output_hdf5_value_cb);
    }
}
void fclaw_clawpatch_output_hdf5 (fclaw_global_t * glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d.hdf5", fclaw_opt->prefix, iframe);

    fclaw_clawpatch_output_hdf5_to_file(glob,basename);
}



