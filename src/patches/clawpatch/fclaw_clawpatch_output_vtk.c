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

#include <fclaw_clawpatch_output_vtk.h>

#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch_vtk_vtable.h>
#include <fclaw_global.h>

/** Round up quotient x / y. */
#define FCLAW_VTK_CEIL(x, y) (((x) + (y) - 1) / (y)) /**< round up the quotient
                                                        x / y; this is a local
                                                        macro used for
                                                        computations regarding
                                                        the VTK patch buffering */

#include <fclaw_options.h>
#include <fclaw_map.h>
#include <sc_io.h>

typedef struct fclaw2d_vtk_state
{
    int dim;
    int patch_children;
    int mx, my, mz;
    int meqn;
    int num_aux_fields;
    int rhs_fields;
    int points_per_patch, cells_per_patch;
    int intsize, ndsize;
    int fits32;
    double time_value;
    char filename[BUFSIZ];
    int64_t global_num_points, global_num_cells;
    int64_t global_num_connectivity;
    int64_t global_num_patches;
    int64_t *field_offsets;
    int64_t *point_offsets;
    int64_t *cell_offsets;
    int64_t *celldata_offsets;
    int64_t offset_end;
    const char *inttype;
    fclaw_vtk_patch_data_t coordinate_cb;
    fclaw_vtk_patch_data_t value_cb;
    fclaw_vtk_patch_data_t aux_cb;
    fclaw_vtk_patch_data_t rhs_cb;
    fclaw_vtk_patch_data_t soln_cb;
    fclaw_vtk_patch_data_t error_cb;
    fclaw_clawpatch_vtk_vtable_t *vtk_vt;
    FILE *file;
#ifdef P4EST_ENABLE_MPIIO
    MPI_File mpifile;
    MPI_Offset mpibegin;
#endif
    char *buf;

    /* The following three struct elements are dedicated to manage buffering
     * patches.
     * We track the number of buffered patches in \b num_buffered_patches and
     * if the buffer contains more than \b patch_threshold many patches, the
     * patches are written to disk and the buffer is flushed.
     * The variable \b num_buffered_patches relates to the number of buffered
     * patches during one fclaw2d_vtk_write_field call.
     * If \b patch_threshold is 0, there are no intermediate writes but all
     * patches during a fclaw2d_vtk_write_field are buffered and then written
     * to the disk.
     */
    sc_io_sink_t *sink; /**< the sink to manage the patches buffer */
    int num_buffered_patches; /**< the number of patches buffered during the
                                   current fclaw2d_vtk_write_field call */
    int patch_threshold; /**< maximal number of patches that are buffered during
                              a call of fclaw2d_vtk_write_field */
}
fclaw2d_vtk_state_t;

static size_t
sizeof_vtk_type(fclaw2d_vtk_state_t *s, fclaw_vtk_entry_type_t type)
{
    switch (type)
    {
        case FCLAW_VTK_UINT8:
            return sizeof(uint8_t);
        case FCLAW_VTK_INT32:
            return sizeof(int32_t);
        case FCLAW_VTK_INT64:
            return sizeof(int64_t);
        case FCLAW_VTK_UINT64:
            return sizeof(uint64_t);
        case FCLAW_VTK_FLOAT32:
            return sizeof(float);
        case FCLAW_VTK_FLOAT64:
            return sizeof(double);
        case FCLAW_VTK_INT32_OR_64:
            /* this is a special entry type that is used for connectivity and
             * offsets; the actual type used depends on whether the number of
             * points/cells fits in 32 bits */
            return s->fits32 ? sizeof(int32_t) : sizeof(int64_t);
        default:
            fclaw_abortf("Invalid VTK entry type %d\n", (int) type);
    }
}

static const char*
vtk_type_to_string(fclaw2d_vtk_state_t *s, fclaw_vtk_entry_type_t type)
{
    switch (type)
    {
        case FCLAW_VTK_UINT8:
            return "UInt8";
        case FCLAW_VTK_INT32:
            return "Int32";
        case FCLAW_VTK_INT64:
            return "Int64";
        case FCLAW_VTK_UINT64:
            return "UInt64";
        case FCLAW_VTK_FLOAT32:
            return "Float32";
        case FCLAW_VTK_FLOAT64:
            return "Float64";
        case FCLAW_VTK_INT32_OR_64:
            return s->fits32 ? "Int32" : "Int64";
        default:
            fclaw_abortf("Invalid VTK entry type %d\n", (int) type);
    }
}

static int
print_dataarray_entry(FILE *file,
                      fclaw2d_vtk_state_t *s,
                      const fclaw_clawpatch_vtk_vtable_entry_t *entry,
                      int64_t offset,
                      const char *format_single,
                      const char *format_multiple)
{
    const char* vtk_type_str = vtk_type_to_string(s, entry->type);
    if (entry->number_of_components == 1)
    {
        return fprintf(file, format_single,
                       vtk_type_str,
                       entry->name,
                       (long long) offset) < 0;
    }
    return fprintf(file, format_multiple,
                   vtk_type_str,
                   entry->name,
                   entry->number_of_components,
                   (long long) offset) < 0;
}

static int
fclaw2d_vtk_write_header (fclaw_domain_t * domain, fclaw2d_vtk_state_t * s)
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
    retval = retval || fprintf (file, "  <FieldData>\n") < 0;
    retval = retval || fprintf (file, "   <DataArray type=\"Int32\" "
                                "Name=\"patch_dimension\" NumberOfComponents=\"3\" NumberOfTuples=\"1\" format=\"ascii\">\n") < 0;
    retval = retval || fprintf (file, "    %d %d %d\n", s->mx, s->my, s->mz) < 0;
    retval = retval || fprintf (file, "   </DataArray>\n") < 0;
    retval = retval || fprintf (file, "   <DataArray type=\"Float64\" "
                                "Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">\n") < 0;
    retval = retval || fprintf (file, "    %.*g\n", 17, s->time_value) < 0;
    retval = retval || fprintf (file, "   </DataArray>\n") < 0;

    const char *field_format_single = "   <DataArray type=\"%s\" Name=\"%s\" NumberOfTuples=\"%lld\" format=\"appended\" offset=\"%lld\">\n"
                                      "   </DataArray>\n";
    const char *field_format_multiple = "   <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" NumberOfTuples=\"%lld\" format=\"appended\" offset=\"%lld\">\n"
                                        "   </DataArray>\n";

    // write field vtable entries
    sc_link_t* curr_entry = s->vtk_vt->field_entries->first;
    int64_t* curr_offset = s->field_offsets;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        const char* vtk_type_str = vtk_type_to_string(s, entry->type);
        long long number_of_tuples = (long long) entry->elements_per_patch * s->global_num_patches;
        if (entry->number_of_components == 1)
        {
            retval = retval || fprintf (file, field_format_single,
                                        vtk_type_str,
                                        entry->name,
                                        (long long) number_of_tuples,
                                        (long long) *curr_offset) < 0;
        }
        else
        {
            retval = retval || fprintf (file, field_format_multiple,
                                        vtk_type_str,
                                        entry->name,
                                        entry->number_of_components,
                                        (long long) number_of_tuples,
                                        (long long) *curr_offset) < 0;
        }
        curr_offset++;
        curr_entry = curr_entry->next;
    }

    retval = retval || fprintf (file, "  </FieldData>\n") < 0;
    retval = retval || fprintf (file, "  <Piece NumberOfPoints=\"%lld\" "
                                "NumberOfCells=\"%lld\">\n",
                                (long long) s->global_num_points,
                                (long long) s->global_num_cells) < 0;
    retval = retval || fprintf (file, "   <Points>\n") < 0;

    const char *format_single = "    <DataArray type=\"%s\" Name=\"%s\" format=\"appended\" offset=\"%lld\">\n"
                                "    </DataArray>\n";
    const char *format_multiple = "    <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\" offset=\"%lld\">\n"
                                  "    </DataArray>\n";

    // write primary point entry
    retval = retval || print_dataarray_entry(file,
                                             s,
                                             &s->vtk_vt->position_entry,
                                             *curr_offset,
                                             format_single,
                                             format_multiple);
    curr_offset++;

    retval = retval || fprintf (file, "   </Points>\n") < 0;
    retval = retval || fprintf (file, "   <Cells>\n") < 0;

    // write primary cell entries
    retval = retval || print_dataarray_entry(file,
                                             s,
                                             &s->vtk_vt->connectivity_entry,
                                             *curr_offset,
                                             format_single,
                                             format_multiple);
    curr_offset++;
    retval = retval || print_dataarray_entry(file,
                                             s,
                                             &s->vtk_vt->offsets_entry,
                                             *curr_offset,
                                             format_single,
                                             format_multiple);
    curr_offset++;
    retval = retval || print_dataarray_entry(file,
                                             s,
                                             &s->vtk_vt->types_entry,
                                             *curr_offset,
                                             format_single,
                                             format_multiple);
    curr_offset++;

    retval = retval || fprintf (file, "   </Cells>\n") < 0;
    retval = retval || fprintf (file, "   <CellData>\n") < 0;

    // write celldata vtable entries
    curr_entry = s->vtk_vt->celldata_entries->first;
    curr_offset = s->celldata_offsets;
    while (curr_entry != NULL)
    {        
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        const char* vtk_type_str = vtk_type_to_string(s, entry->type);
        if (entry->number_of_components == 1)
        {            retval = retval || fprintf (file, format_single,
                                        vtk_type_str,
                                        entry->name,
                                        (long long) *curr_offset) < 0;
        }
        else        
        {
            retval = retval || fprintf (file, format_multiple,
                                        vtk_type_str,
                                        entry->name,
                                        entry->number_of_components,
                                        (long long) *curr_offset) < 0;
        }
        curr_offset++;
        curr_entry = curr_entry->next;
    }

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

/** This function adds to a buffer and writes to file if a threshold is exceeded.
 *
 * This function assumes that the buffer consists of patch data of the size
 * \b psize_field. This means it must hold
 * \b s->sink->buffer_bytes % \b psize_field == 0.
 *
 * \param [in,out]  s               The VTK state.
 * \param [in]      psize_field     The number of bytes, which are intended to
 *                                  add to the buffer.
 * \param [in]      threshold       If the number of patches, which are already
 *                                  stored in the buffer + 1 is greater than
 *                                  \b threshold, the buffer is flushed to disk.
 *                                  Then the new patch is added to the buffer.
 *                                  0 means that there is no threshold and
 *                                  the data just added to \b s->sink.
 */
static void
add_to_buffer (fclaw2d_vtk_state_t * s, int64_t psize_field)
{
    FCLAW_ASSERT (s != NULL);

    int mpiret;
#ifdef P4EST_ENABLE_MPIIO
    MPI_Status mpistatus;
#else
    size_t retvalz;
#endif
    FCLAW_ASSERT (s->patch_threshold >= 0);
    FCLAW_ASSERT (((int) s->sink->buffer_bytes) % psize_field == 0);

    if (s->patch_threshold != 0
        && (s->num_buffered_patches + 1 > s->patch_threshold))
    {
        FCLAW_ASSERT ((int) s->sink->buffer_bytes >= psize_field);
        FCLAW_ASSERT (((int) s->sink->buffer_bytes) / psize_field ==
                      s->num_buffered_patches);
        /* buffer threshold exceeded */
        /* write the current buffer to disk */
#ifdef P4EST_ENABLE_MPIIO
        if (s->num_buffered_patches > 0)
        {
            mpiret =
                MPI_File_write_all (s->mpifile, s->sink->buffer->array,
                                    (int) s->sink->buffer_bytes, MPI_BYTE,
                                    &mpistatus);
            SC_CHECK_MPI (mpiret);
        }
#else
        retvalz =
            fwrite (s->sink->buffer->array, s->sink->buffer_bytes, 1,
                    s->file);
        SC_CHECK_ABORT (retvalz == 1, "VTK file write failed");
#endif
        /* reset the sink */
        mpiret = sc_io_sink_complete (s->sink, NULL, NULL);
        SC_CHECK_ABORT (mpiret == 0, "VTK buffer completion failed");
        /* reset the buffer */
        sc_array_reset (s->sink->buffer);
        s->sink->buffer_bytes = 0;
        s->num_buffered_patches = 0;
    }

    mpiret = sc_io_sink_write (s->sink, s->buf, (size_t) psize_field);
    SC_CHECK_ABORT (mpiret == 0, "VTK buffer write failed");
    ++s->num_buffered_patches;
}

typedef struct write_field_iter_user
{
    fclaw2d_vtk_state_t *s; /**< the VTK state */
    fclaw_vtk_cb_context_t ctx; /**< the VTK callback context */
    fclaw_clawpatch_vtk_vtable_entry_t* entry; /**< the vtable entry for the field */
} write_field_iter_user_t;

static void
write_entry_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    write_field_iter_user_t *iter = (write_field_iter_user_t *) g->user;

    iter->entry->callback (g->glob, patch, blockno, patchno, &iter->ctx, iter->s->buf);

    int64_t psize = iter->entry->elements_per_patch * iter->entry->number_of_components * sizeof_vtk_type(iter->s, iter->entry->type);
    add_to_buffer (iter->s, psize);
}

static void
fclaw2d_vtk_write_field (fclaw_global_t * glob, fclaw2d_vtk_state_t * s,
                         int64_t offset_field, fclaw_clawpatch_vtk_vtable_entry_t* entry)
{
    int64_t psize_field = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    if(psize_field == 0)
    {
        /* some fields are optional */
        return;
    }

    fclaw_domain_t *domain = glob->domain;

    int64_t bcount;
    sc_array_t buffer;
#ifndef P4EST_ENABLE_MPIIO
    size_t retvalz;
#else
    int i;
    int max_num_writes, local_num_writes;
    int mpiret;
    MPI_Offset mpipos;
#ifdef FCLAW_ENABLE_DEBUG
    MPI_Offset mpinew;
#endif
    MPI_Status mpistatus;
#endif

    s->buf = P4EST_ALLOC (char, psize_field);
    /* The buffer array will be resized when required. */
    sc_array_init_size (&buffer, 1, psize_field);
    s->sink = sc_io_sink_new (SC_IO_TYPE_BUFFER, SC_IO_MODE_APPEND,
                              SC_IO_ENCODE_NONE, &buffer);
    FCLAW_ASSERT (s->num_buffered_patches == 0);
    s->num_buffered_patches = 0;
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
        mpiret =
            MPI_File_write (s->mpifile, &bcount, 1, MPI_LONG, &mpistatus);
        SC_CHECK_MPI (mpiret);
#endif
    }
    write_field_iter_user_t iter;
    iter.s = s;
    iter.ctx.fits32 = s->fits32;
    iter.entry = entry;
    fclaw_global_iterate_patches (glob, write_entry_cb, &iter);

#ifdef P4EST_ENABLE_MPIIO
    if (s->num_buffered_patches > 0 || s->patch_threshold == 0)
    {
        /* write the remaining buffered bytes */
        mpiret =
            MPI_File_write_all (s->mpifile, buffer.array,
                                (int) s->sink->buffer_bytes, MPI_BYTE,
                                &mpistatus);
        SC_CHECK_MPI (mpiret);
    }

    if (s->patch_threshold != 0)
    {
        /* Ensure that the collective function MPI_File_write_all is called
         * equally often on each rank.
         */
        max_num_writes =
            FCLAW_VTK_CEIL (glob->domain->local_max_patches,
                            s->patch_threshold);
        local_num_writes =
            FCLAW_VTK_CEIL (glob->domain->local_num_patches,
                            s->patch_threshold);
        FCLAW_ASSERT (max_num_writes - local_num_writes >= 0);

        for (i = 0; i < max_num_writes - local_num_writes; ++i)
        {
            /* This rank has less patches than the rank that holds the maximal
             * number of patches. We compensate this by empty write calls to avoid
             * a deadlock.
             */
            mpiret =
                MPI_File_write_all (s->mpifile, buffer.array, 0, MPI_BYTE,
                                    &mpistatus);
            SC_CHECK_MPI (mpiret);
        }
    }
#else
    retvalz = fwrite (buffer.array, s->sink->buffer_bytes, 1, s->file);
    SC_CHECK_ABORT (retvalz == 1, "VTK file write failed");
#endif

    /* free buffers */
    sc_array_reset (&buffer);
    sc_io_sink_destroy (s->sink);
    P4EST_FREE (s->buf);
    s->num_buffered_patches = 0;

#ifdef P4EST_ENABLE_MPIIO
#ifdef FCLAW_ENABLE_DEBUG
    mpiret = MPI_File_get_position (s->mpifile, &mpinew);
    SC_CHECK_MPI (mpiret);
    FCLAW_ASSERT (mpinew - mpipos ==
                  (domain->mpirank == 0 ? s->ndsize : 0) +
                  domain->local_num_patches * psize_field);
    FCLAW_ASSERT (domain->mpirank < domain->mpisize - 1 ||
                  mpinew - s->mpibegin ==
                  offset_field + s->ndsize +
                  psize_field * domain->global_num_patches);
#endif
#endif
}

static void
fclaw2d_vtk_write_data (fclaw_global_t * glob, fclaw2d_vtk_state_t * s)
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

    // write field vtable entries
    sc_link_t* curr_entry = s->vtk_vt->field_entries->first;
    int64_t* curr_offset = s->field_offsets;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        fclaw2d_vtk_write_field (glob, s, *curr_offset, entry);
        curr_offset++;
        curr_entry = curr_entry->next;
    }

    // write primary point entry
    curr_offset = s->point_offsets;
    fclaw2d_vtk_write_field(glob, s, *curr_offset, &s->vtk_vt->position_entry);
    curr_offset++;

    // write primary cell entries
    curr_offset = s->cell_offsets;
    fclaw2d_vtk_write_field(glob, s, *curr_offset, &s->vtk_vt->connectivity_entry);
    curr_offset++;
    fclaw2d_vtk_write_field(glob, s, *curr_offset, &s->vtk_vt->offsets_entry);
    curr_offset++;
    fclaw2d_vtk_write_field(glob, s, *curr_offset, &s->vtk_vt->types_entry);
    curr_offset++;

    // write celldata vtable entries
    curr_entry = s->vtk_vt->celldata_entries->first;
    curr_offset = s->celldata_offsets;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        fclaw2d_vtk_write_field (glob, s, *curr_offset, entry);
        curr_offset++;
        curr_entry = curr_entry->next;
    }

#ifdef P4EST_ENABLE_MPIIO
    /* collectively close the file */
    mpiret = MPI_File_close (&s->mpifile);
    SC_CHECK_MPI (mpiret);
#endif
}

static int
fclaw2d_vtk_write_footer (fclaw_domain_t * domain, fclaw2d_vtk_state_t * s)
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

static int
fclaw_vtk_write_file (int dim, fclaw_global_t * glob, const char *basename,
                      int mx, int my, int mz,
                      int meqn, int num_aux_fields, int rhs_fields,
                      double vtkspace, int vtkwrite,
                      fclaw_vtk_patch_data_t coordinate_cb,
                      fclaw_vtk_patch_data_t value_cb,
                      fclaw_vtk_patch_data_t aux_cb,
                      fclaw_vtk_patch_data_t rhs_cb,
                      fclaw_vtk_patch_data_t soln_cb,
                      fclaw_vtk_patch_data_t error_cb,
                      int patch_threshold)
{
    fclaw_domain_t *domain = glob->domain;

    FCLAW_ASSERT (patch_threshold >= 0);

    int retval, gretval;
    int mpiret;
    fclaw2d_vtk_state_t ps, *s = &ps;

    /* set up VTK internal information */
    s->dim = dim;
    s->patch_children = (dim == 2) ? 4 : 8;
    s->mx = mx;
    s->my = my;
    if(dim == 3)
    {
        s->mz = mz;
    }
    else
    {
        s->mz = 1;
    }
    s->meqn = meqn;
    s->rhs_fields = rhs_fields;
    s->num_aux_fields = num_aux_fields;
    s->points_per_patch = (mx + 1) * (my + 1);
    s->cells_per_patch = mx * my;
    if(dim == 3)
    {
        s->points_per_patch *= (mz + 1);
        s->cells_per_patch *= mz;
    }
    snprintf (s->filename, BUFSIZ, "%s.vtu", basename);
    s->global_num_points = s->points_per_patch * domain->global_num_patches;
    s->global_num_cells = s->cells_per_patch * domain->global_num_patches;
    s->global_num_connectivity = s->patch_children * (s->global_num_cells + 1);
    s->global_num_patches = domain->global_num_patches;
    s->fits32 = s->global_num_points <= INT32_MAX
        && s->global_num_connectivity <= INT32_MAX;
    s->inttype = s->fits32 ? "Int32" : "Int64";
    s->intsize = s->fits32 ? sizeof (int32_t) : sizeof (int64_t);
    s->ndsize = 8;   /* uint64 */
    s->time_value = glob->curr_time;
    s->coordinate_cb = coordinate_cb;
    s->value_cb = value_cb;
    s->aux_cb = aux_cb;
    s->rhs_cb = rhs_cb;
    s->soln_cb = soln_cb;
    s->error_cb = error_cb;

    s->vtk_vt = fclaw_clawpatch_vtk_vtable(glob);

    /* compute offsets in bytes after beginning of appended data section */
    int64_t curr_offset = 0;

    int num_field_entries = s->vtk_vt->field_entries->elem_count;
    s->field_offsets = FCLAW_ALLOC (int64_t, num_field_entries);

    sc_link_t* curr_entry = s->vtk_vt->field_entries->first;
    int64_t* curr_entry_offset = s->field_offsets;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        *curr_entry_offset = curr_offset;
        size_t psize = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
        curr_offset += s->ndsize + psize * s->global_num_patches;
        curr_entry_offset++;
        curr_entry = curr_entry->next;
    }


    int num_point_entries = 1;
    s->point_offsets = FCLAW_ALLOC (int64_t, num_point_entries);

    curr_entry_offset = s->point_offsets;

    fclaw_clawpatch_vtk_vtable_entry_t* entry = &s->vtk_vt->position_entry;
    *curr_entry_offset = curr_offset;
    size_t psize = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + psize * s->global_num_patches;
    curr_entry_offset++;

    int num_cell_entries = 3;
    s->cell_offsets = FCLAW_ALLOC (int64_t, num_cell_entries);

    curr_entry_offset = s->cell_offsets;

    entry = &s->vtk_vt->connectivity_entry;
    *curr_entry_offset = curr_offset;
    psize = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + psize * s->global_num_patches;
    curr_entry_offset++;

    entry = &s->vtk_vt->offsets_entry;
    *curr_entry_offset = curr_offset;
    psize = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + psize * s->global_num_patches;
    curr_entry_offset++;

    entry = &s->vtk_vt->types_entry;
    *curr_entry_offset = curr_offset;
    psize = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + psize * s->global_num_patches;
    curr_entry_offset++;


    int num_celldata_entries = s->vtk_vt->celldata_entries->elem_count;
    s->celldata_offsets = FCLAW_ALLOC (int64_t, num_celldata_entries);

    curr_entry = s->vtk_vt->celldata_entries->first;
    curr_entry_offset = s->celldata_offsets;
    while (curr_entry != NULL)    
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        *curr_entry_offset = curr_offset;
        size_t psize = entry->elements_per_patch * entry->number_of_components * sizeof_vtk_type(s, entry->type);
        curr_offset += s->ndsize + psize * s->global_num_patches;
        curr_entry_offset++;
        curr_entry = curr_entry->next;
    }

    s->offset_end = curr_offset;

    s->buf = NULL;
    s->sink = NULL;
    /* See the documentation of fclaw2d_vtk_state_t for further information. */
    s->patch_threshold = patch_threshold;
    s->num_buffered_patches = 0;

    /* write header meta data and check for error */
    retval = 0;
    if (domain->mpirank == 0)
    {
        retval = fclaw2d_vtk_write_header (glob->domain, s);
    }
    mpiret = sc_MPI_Allreduce (&retval, &gretval, 1, sc_MPI_INT, sc_MPI_MIN,
                               domain->mpicomm);
    SC_CHECK_MPI (mpiret);
    if (gretval < 0)
    {
        return -1;
    }

    /* write mesh and numerical data using MPI I/O */
    fclaw2d_vtk_write_data (glob, s);

    /* write footer information and check for error */
    retval = 0;
    if (domain->mpirank == 0)
    {
        retval = fclaw2d_vtk_write_footer (domain, s);
    }

    // free offset arrays
    FCLAW_FREE (s->field_offsets);
    FCLAW_FREE (s->point_offsets);
    FCLAW_FREE (s->cell_offsets);
    FCLAW_FREE (s->celldata_offsets);

    mpiret = sc_MPI_Allreduce (&retval, &gretval, 1, sc_MPI_INT, sc_MPI_MIN,
                               domain->mpicomm);
    SC_CHECK_MPI (mpiret);
    if (gretval < 0)
    {
        return -1;
    }

    return 0;
}

int
fclaw_vtk_write_2d_file (fclaw_global_t * glob, const char *basename,
                        int mx, int my,
                        int meqn,
                        double vtkspace, int vtkwrite,
                        fclaw_vtk_patch_data_t coordinate_cb,
                        fclaw_vtk_patch_data_t value_cb,
                        int patch_threshold)
{
    return fclaw_vtk_write_file(2,glob,basename,mx,my,0,meqn,0,0,vtkspace,vtkwrite,
                                coordinate_cb,value_cb,NULL,NULL,NULL,NULL, patch_threshold);
}

int
fclaw_vtk_write_3d_file (fclaw_global_t * glob, const char *basename,
                        int mx, int my, int mz,
                        int meqn,
                        double vtkspace, int vtkwrite,
                        fclaw_vtk_patch_data_t coordinate_cb,
                        fclaw_vtk_patch_data_t value_cb,
                        int patch_threshold)
{
    return fclaw_vtk_write_file(3,glob,basename,mx,my,mz,meqn,0,0,vtkspace,vtkwrite,
                                coordinate_cb,value_cb,NULL,NULL,NULL,NULL, patch_threshold);
}

static void
fclaw2d_output_vtk_coordinate_cb (fclaw_global_t * glob,
                                  fclaw_patch_t * patch,
                                  int blockno, int patchno,
                                  char *a)
{
    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fclaw_map_context_t* cont = fclaw_map_get(glob);

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
fclaw3d_output_vtk_coordinate_cb (fclaw_global_t * glob,
                                  fclaw_patch_t * patch,
                                  int blockno, int patchno,
                                  char *a)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;
    fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                &xlower,&ylower,&zlower, &dx,&dy, &dz);

    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    fclaw_map_context_t* cont = fclaw_map_get(glob);
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


/* typdef for data acess functions q,rhs,etc */
typedef void (*patch_data_access_t)(struct fclaw_global *glob, struct fclaw_patch *patch, double **q, int *meqn);

static void
pack_data (fclaw_global_t * glob,
           fclaw_patch_t * patch,
           int blockno, int patchno,
           patch_data_access_t access,
           char *a)
{

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
}

static void
fclaw_output_vtk_value_cb (fclaw_global_t * glob,
                           fclaw_patch_t * patch,
                           int blockno, int patchno,
                           char *a)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_soln_data, a);
}

static void
fclaw_output_vtk_aux_cb (fclaw_global_t * glob,
                         fclaw_patch_t * patch,
                         int blockno, int patchno,
                         char *a)
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
        float *f = (float *) a;
        int i, j, k;
        for (j = 0; j < my; ++j)
        {
            for (i = 0; i < mx; ++i)
            {
                for (k = 0; k < maux; ++k)
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
        float *f = (float *) a;
        int i, j, k, eqn;
        for (k = 0; k < mz; ++k)
        {
            for (j = 0; j < my; ++j)
            {
                for (i = 0; i < mx; ++i)
                {
                    for (eqn = 0; eqn < maux; ++eqn)
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
fclaw_output_vtk_rhs_cb (fclaw_global_t * glob,
                         fclaw_patch_t * patch,
                         int blockno, int patchno,
                         char *a)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_rhs_data, a);
}

static void
fclaw_output_vtk_soln_cb (fclaw_global_t * glob,
                          fclaw_patch_t * patch,
                          int blockno, int patchno,
                          char *a)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_elliptic_soln_data, a);
}

static void
fclaw_output_vtk_error_cb (fclaw_global_t * glob,
                           fclaw_patch_t * patch,
                           int blockno, int patchno,
                           char *a)
{
    pack_data(glob, patch, blockno, patchno, &fclaw_clawpatch_elliptic_error_data, a);
}

/*  --------------------------------------------------------------------------
    Used for debugging
    ------------------------------------------------------------------------- */
#if 0
static void
fclaw2d_output_write_vtk_debug (fclaw_global_t * glob, const char *basename)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    (void) fclaw2d_vtk_write_file (glob, basename,
                                   clawpatch_opt->mx, clawpatch_opt->my,
#if PATCH_DIM == 3
                                   clawpatch_opt->mz,
#endif
                                   clawpatch_opt->meqn,
                                   fclaw_opt->vtkspace, 0,
                                   fclaw2d_output_vtk_coordinate_cb,
                                   fclaw2d_output_vtk_value_cb);
}
#endif


/*  ---------------------------------------------------------------------------
    Public interface
    --------------------------------------------------------------------------- */

void fclaw_clawpatch_output_vtk_to_file (fclaw_global_t * glob, const char* filename)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int num_aux_fields = 0;
    for(int i = 0; i < clawpatch_opt->maux; i++)
    {
        if (clawpatch_opt->vtk_aux_out[i] > 0)
        {
            num_aux_fields++;
        }
    }

    if(clawpatch_opt->patch_dim == 2)
    {
        fclaw_vtk_write_file (2, glob, filename,
                              clawpatch_opt->mx, 
                              clawpatch_opt->my,
                              0,
                              clawpatch_opt->meqn,
                              num_aux_fields,
                              clawpatch_opt->rhs_fields,
                              fclaw_opt->vtkspace, 0,
                              fclaw2d_output_vtk_coordinate_cb,
                              fclaw_output_vtk_value_cb,
                              fclaw_output_vtk_aux_cb,
                              fclaw_output_vtk_rhs_cb,
                              fclaw_output_vtk_soln_cb,
                              fclaw_output_vtk_error_cb,
                              clawpatch_opt->vtk_patch_threshold);
    }
    else 
    {
        fclaw_vtk_write_file (3, glob, filename,
                              clawpatch_opt->mx, 
                              clawpatch_opt->my, 
                              clawpatch_opt->mz,
                              clawpatch_opt->meqn,
                              num_aux_fields,
                              clawpatch_opt->rhs_fields,
                              fclaw_opt->vtkspace, 0,
                              fclaw3d_output_vtk_coordinate_cb,
                              fclaw_output_vtk_value_cb,
                              fclaw_output_vtk_aux_cb,
                              fclaw_output_vtk_rhs_cb,
                              fclaw_output_vtk_soln_cb,
                              fclaw_output_vtk_error_cb,
                              clawpatch_opt->vtk_patch_threshold);
    }
}
void fclaw_clawpatch_output_vtk (fclaw_global_t * glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d", fclaw_opt->prefix, iframe);

    fclaw_clawpatch_output_vtk_to_file(glob,basename);
}



