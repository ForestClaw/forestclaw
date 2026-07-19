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

#include <fclaw_config.h>
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
    int64_t *field_num_elements;
    int64_t *point_num_elements;
    int64_t *cell_num_elements;
    int64_t *celldata_num_elements;
    int64_t offset_end;
    const char *inttype;
    fclaw_vtk_patch_data_t coordinate_cb;
    fclaw_vtk_patch_data_t value_cb;
    fclaw_clawpatch_vtk_vtable_t *vtk_vt;
    FILE *file;
#ifdef P4EST_ENABLE_MPIIO
    MPI_File mpifile;
    MPI_Offset mpibegin;
#endif
    char *buf;
    size_t buf_capacity;

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

typedef struct count_elements_user
{
    fclaw2d_vtk_state_t *s;
    fclaw_vtk_cb_context_t ctx;
    fclaw_clawpatch_vtk_vtable_entry_t *entry;
    unsigned long long num_elements;
} count_elements_user_t;

static void
count_elements_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                   int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    count_elements_user_t *iter = (count_elements_user_t *) g->user;

    iter->num_elements += iter->entry->elements_in_patch (g->glob, patch, blockno, patchno, &iter->ctx);
}

static int64_t
get_global_num_elements (fclaw_global_t *glob, fclaw2d_vtk_state_t *s, fclaw_clawpatch_vtk_vtable_entry_t *entry)
{
    fclaw_domain_t *domain = glob->domain;

    count_elements_user_t count_user;
    count_user.s = s;
    count_user.ctx.fits32 = s->fits32;
    count_user.ctx.offsets_include_zero = 0;
    count_user.entry = entry;
    count_user.num_elements = 0;

    fclaw_global_iterate_patches (glob, count_elements_cb, &count_user);

    unsigned long long local_num_elements = count_user.num_elements;
    unsigned long long global_num_elements = 0;
    int mpiret = sc_MPI_Allreduce (&local_num_elements, &global_num_elements, 1,
                                   sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM, domain->mpicomm);
    SC_CHECK_MPI (mpiret);

    return (int64_t) global_num_elements;
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
    int i = 0;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        const char* vtk_type_str = vtk_type_to_string(s, entry->type);
        long long number_of_tuples = (long long) s->field_num_elements[i];
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
        i++;
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
    curr_offset = s->point_offsets;
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
    curr_offset = s->cell_offsets;
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

    if (s->patch_threshold != 0
        && (s->num_buffered_patches + 1 > s->patch_threshold))
    {
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

    size_t num_elements = iter->entry->elements_in_patch (g->glob, patch, blockno, patchno, &iter->ctx);
    int64_t psize = num_elements * iter->entry->number_of_components * sizeof_vtk_type(iter->s, iter->entry->type);

    if (psize > 0)
    {
        if ((size_t) psize > iter->s->buf_capacity)
        {
            if (iter->s->buf != NULL)
            {
                P4EST_FREE (iter->s->buf);
            }
            iter->s->buf = P4EST_ALLOC (char, psize);
            iter->s->buf_capacity = psize;
        }

        if (iter->entry == &iter->s->vtk_vt->position_entry && iter->s->coordinate_cb != NULL)
        {
            iter->s->coordinate_cb (g->glob, patch, blockno, patchno, iter->s->buf);
        }
        else if (strcmp(iter->entry->name, "meqn") == 0 && iter->s->value_cb != NULL)
        {
            iter->s->value_cb (g->glob, patch, blockno, patchno, iter->s->buf);
        }
        else
        {
            iter->entry->callback (g->glob, patch, blockno, patchno, &iter->ctx, iter->s->buf);
        }
        add_to_buffer (iter->s, psize);
    }
}

typedef struct count_bytes_user
{
    fclaw2d_vtk_state_t *s;
    fclaw_vtk_cb_context_t ctx;
    fclaw_clawpatch_vtk_vtable_entry_t *entry;
    int64_t local_bytes;
} count_bytes_user_t;

static void
count_bytes_cb (fclaw_domain_t * domain, fclaw_patch_t * patch,
                int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t *g = (fclaw_global_iterate_t*) user;
    count_bytes_user_t *iter = (count_bytes_user_t *) g->user;

    size_t num_elements = iter->entry->elements_in_patch (g->glob, patch, blockno, patchno, &iter->ctx);
    int64_t psize = num_elements * iter->entry->number_of_components * sizeof_vtk_type(iter->s, iter->entry->type);
    iter->local_bytes += psize;
}

static void
fclaw2d_vtk_write_field (fclaw_global_t * glob, fclaw2d_vtk_state_t * s,
                         int64_t offset_field, fclaw_clawpatch_vtk_vtable_entry_t* entry)
{
    int64_t global_num_elements = get_global_num_elements (glob, s, entry);
    if (global_num_elements == 0)
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

    count_bytes_user_t bytes_user;
    bytes_user.s = s;
    bytes_user.ctx.fits32 = s->fits32;
    bytes_user.ctx.offsets_include_zero = 0;
    bytes_user.entry = entry;
    bytes_user.local_bytes = 0;

    fclaw_global_iterate_patches (glob, count_bytes_cb, &bytes_user);

    unsigned long long local_bytes = bytes_user.local_bytes;
    unsigned long long prefix_bytes = 0;
    unsigned long long global_bytes = 0;
#ifdef FCLAW_ENABLE_MPI
    int mpiret_allreduce = sc_MPI_Scan (&local_bytes, &prefix_bytes, 1, sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM, domain->mpicomm);
    SC_CHECK_MPI (mpiret_allreduce);
    mpiret_allreduce = sc_MPI_Allreduce (&local_bytes, &global_bytes, 1, sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM, domain->mpicomm);
    SC_CHECK_MPI (mpiret_allreduce);
#else
    prefix_bytes = local_bytes;
    global_bytes = local_bytes;
#endif

    unsigned long long bytes_before = prefix_bytes - local_bytes;

    s->buf = NULL;
    s->buf_capacity = 0;

    /* The buffer array will be resized when required. */
    sc_array_init_size (&buffer, 1, 0);
    s->sink = sc_io_sink_new (SC_IO_TYPE_BUFFER, SC_IO_MODE_APPEND,
                              SC_IO_ENCODE_NONE, &buffer);
    FCLAW_ASSERT (s->num_buffered_patches == 0);
    s->num_buffered_patches = 0;
#ifdef P4EST_ENABLE_MPIIO
    mpipos = s->mpibegin + offset_field;
    if (domain->mpirank > 0)
    {
        /* account for byte count */
        mpipos += s->ndsize + bytes_before;
    }
    mpiret = MPI_File_seek (s->mpifile, mpipos, MPI_SEEK_SET);
    SC_CHECK_MPI (mpiret);
#endif
    if (domain->mpirank == 0)
    {
        /* write byte count */
        bcount = (int64_t) global_bytes;
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
    iter.ctx.offsets_include_zero = 0;
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
    if (s->buf != NULL)
    {
        P4EST_FREE (s->buf);
        s->buf = NULL;
    }
    s->buf_capacity = 0;
    s->num_buffered_patches = 0;

#ifdef P4EST_ENABLE_MPIIO
#ifdef FCLAW_ENABLE_DEBUG
    mpiret = MPI_File_get_position (s->mpifile, &mpinew);
    SC_CHECK_MPI (mpiret);
    FCLAW_ASSERT (mpinew - mpipos ==
                  (domain->mpirank == 0 ? s->ndsize : 0) +
                  local_bytes);
    FCLAW_ASSERT (domain->mpirank < domain->mpisize - 1 ||
                  mpinew - s->mpibegin ==
                  offset_field + s->ndsize +
                  global_bytes);
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
    s->points_per_patch = 0; /* not used */
    s->cells_per_patch = 0;  /* not used */
    snprintf (s->filename, BUFSIZ, "%s.vtu", basename);

    s->vtk_vt = fclaw_clawpatch_vtk_vtable(glob);

    s->fits32 = 1;
    s->global_num_points = get_global_num_elements (glob, s, &s->vtk_vt->position_entry);
    s->global_num_cells = get_global_num_elements (glob, s, &s->vtk_vt->types_entry);
    s->global_num_connectivity = get_global_num_elements (glob, s, &s->vtk_vt->connectivity_entry);

    s->global_num_patches = domain->global_num_patches;
    s->fits32 = s->global_num_points <= INT32_MAX
        && s->global_num_connectivity <= INT32_MAX;
    s->inttype = s->fits32 ? "Int32" : "Int64";
    s->intsize = s->fits32 ? sizeof (int32_t) : sizeof (int64_t);
    s->ndsize = 8;   /* uint64 */
    s->time_value = glob->curr_time;
    s->coordinate_cb = coordinate_cb;
    s->value_cb = value_cb;

    /* compute offsets in bytes after beginning of appended data section */
    int64_t curr_offset = 0;

    int num_field_entries = s->vtk_vt->field_entries->elem_count;
    s->field_offsets = FCLAW_ALLOC (int64_t, num_field_entries);
    s->field_num_elements = FCLAW_ALLOC (int64_t, num_field_entries);

    sc_link_t* curr_entry = s->vtk_vt->field_entries->first;
    int64_t* curr_entry_offset = s->field_offsets;
    int i = 0;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        *curr_entry_offset = curr_offset;
        int64_t global_num_elements = get_global_num_elements (glob, s, entry);
        s->field_num_elements[i] = global_num_elements;
        int64_t total_bytes = global_num_elements * entry->number_of_components * sizeof_vtk_type(s, entry->type);
        curr_offset += s->ndsize + total_bytes;
        curr_entry_offset++;
        curr_entry = curr_entry->next;
        i++;
    }


    int num_point_entries = 1;
    s->point_offsets = FCLAW_ALLOC (int64_t, num_point_entries);
    s->point_num_elements = FCLAW_ALLOC (int64_t, num_point_entries);

    curr_entry_offset = s->point_offsets;

    fclaw_clawpatch_vtk_vtable_entry_t* entry = &s->vtk_vt->position_entry;
    *curr_entry_offset = curr_offset;
    int64_t global_num_elements = get_global_num_elements (glob, s, entry);
    s->point_num_elements[0] = global_num_elements;
    int64_t total_bytes = global_num_elements * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + total_bytes;
    curr_entry_offset++;

    int num_cell_entries = 3;
    s->cell_offsets = FCLAW_ALLOC (int64_t, num_cell_entries);
    s->cell_num_elements = FCLAW_ALLOC (int64_t, num_cell_entries);

    curr_entry_offset = s->cell_offsets;

    entry = &s->vtk_vt->connectivity_entry;
    *curr_entry_offset = curr_offset;
    global_num_elements = get_global_num_elements (glob, s, entry);
    s->cell_num_elements[0] = global_num_elements;
    total_bytes = global_num_elements * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + total_bytes;
    curr_entry_offset++;

    entry = &s->vtk_vt->offsets_entry;
    *curr_entry_offset = curr_offset;
    global_num_elements = get_global_num_elements (glob, s, entry);
    s->cell_num_elements[1] = global_num_elements;
    total_bytes = global_num_elements * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + total_bytes;
    curr_entry_offset++;

    entry = &s->vtk_vt->types_entry;
    *curr_entry_offset = curr_offset;
    global_num_elements = get_global_num_elements (glob, s, entry);
    s->cell_num_elements[2] = global_num_elements;
    total_bytes = global_num_elements * entry->number_of_components * sizeof_vtk_type(s, entry->type);
    curr_offset += s->ndsize + total_bytes;
    curr_entry_offset++;


    int num_celldata_entries = s->vtk_vt->celldata_entries->elem_count;
    s->celldata_offsets = FCLAW_ALLOC (int64_t, num_celldata_entries);
    s->celldata_num_elements = FCLAW_ALLOC (int64_t, num_celldata_entries);

    curr_entry = s->vtk_vt->celldata_entries->first;
    curr_entry_offset = s->celldata_offsets;
    i = 0;
    while (curr_entry != NULL)    
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        *curr_entry_offset = curr_offset;
        global_num_elements = get_global_num_elements (glob, s, entry);
        s->celldata_num_elements[i] = global_num_elements;
        total_bytes = global_num_elements * entry->number_of_components * sizeof_vtk_type(s, entry->type);
        curr_offset += s->ndsize + total_bytes;
        curr_entry_offset++;
        curr_entry = curr_entry->next;
        i++;
    }

    s->offset_end = curr_offset;

    s->buf = NULL;
    s->buf_capacity = 0;
    s->sink = NULL;
    /* See the documentation of fclaw2d_vtk_state_t for further information. */
    s->patch_threshold = patch_threshold;
    s->num_buffered_patches = 0;

    /* write header meta data and check for error */
    retval = 0;
    if (domain->mpirank == 0)
    {
        retval = fclaw2d_vtk_write_header (domain, s);
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

    // free offset and count arrays
    FCLAW_FREE (s->field_offsets);
    FCLAW_FREE (s->point_offsets);
    FCLAW_FREE (s->cell_offsets);
    FCLAW_FREE (s->celldata_offsets);
    FCLAW_FREE (s->field_num_elements);
    FCLAW_FREE (s->point_num_elements);
    FCLAW_FREE (s->cell_num_elements);
    FCLAW_FREE (s->celldata_num_elements);

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
                                coordinate_cb,value_cb, patch_threshold);
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
                                coordinate_cb,value_cb, patch_threshold);
}

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
                              NULL,
                              NULL,
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
                              NULL,
                              NULL,
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



