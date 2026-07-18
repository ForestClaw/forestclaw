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
#include <fclaw_clawpatch_vtk_vtable.h>

#include <fclaw_global.h>
#include <fclaw_filesystem.h>

#include <fclaw_options.h>
#include <hdf5.h>
#include <string.h>

/*----------------------------------------------------------------------
    Utility functions
  ----------------------------------------------------------------------*/

/**
 * Retrieves the block number and patch number for a given index.
 *
 * @param glob The global context.
 * @param i The patch index.
 * @param blockno Pointer to store the block number.
 * @param patchno Pointer to store the patch number.
 */
static void
get_patch_blockno(fclaw_global_t *glob, int i, int *blockno, int *patchno)
{
    *blockno = 0;
    while(*blockno < (glob->domain->num_blocks-1) && glob->domain->blocks[*blockno+1].num_patches_before <= i)
    {
        (*blockno)++;
    }
    *patchno = i - glob->domain->blocks[*blockno].num_patches_before;
}

/*----------------------------------------------------------------------
    HDF funcitons
  ----------------------------------------------------------------------*/

/**
 * Sets a numerical attribute for an HDF5 object.
 *
 * @param loc_id The location identifier of the object.
 * @param obj_name The name of the object.
 * @param attr_name The name of the attribute.
 * @param size The size of the attribute data.
 * @param tid The datatype the attribute.
 * @param data The attribute data.
 */
static
void set_attribute_numerical(hid_t loc_id, 
                             const char *obj_name, 
                             const char *attr_name, 
                             hsize_t size,
                             hid_t tid, 
                             const void *data)
{
    herr_t status = 0;

    /* Open the object */
    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);

    /* Create the data space for the attribute. */
    hid_t sid = H5Screate_simple(1, &size, NULL);

    /* Delete the attribute if it already exists */
    if(H5Aexists(obj_id, attr_name))
        status |= H5Adelete(obj_id, attr_name);

    /* Create the attribute. */
    hid_t attr_id = H5Acreate2(obj_id, attr_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT);

    /* Write the attribute data. */
    status |= H5Awrite(attr_id, tid, data);

    status |= H5Aclose(attr_id);
    status |= H5Sclose(sid);
    status |= H5Oclose(obj_id);

    if(status != 0 || obj_id < 0 || sid < 0 || attr_id < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in set_attribute_numerical\n");
    }
}

/**
 * Creates a single-value dataset in an HDF5 file.
 *
 * This is used for vtu numcells, numconnectivity, etc
 *
 * @param mpirank The MPI rank of the current process.
 * @param loc_id The identifier of the location where the dataset will be created.
 * @param dset_name The name of the dataset.
 * @param rank The rank of the dataset.
 * @param dims An array of dimensions specifying the size of the dataset.
 * @param tid The datatype of the dataset.
 * @param data A pointer to the data to be written to the dataset.
 */
static void make_single_value_dataset_numerical(int mpirank,
                                                hid_t loc_id,
                                                const char *dset_name,
                                                int rank,
                                                const hsize_t *dims,
                                                hid_t tid,
                                                const void *data)
{
    herr_t status = 0;

    /* Create the data space for the dataset. */
    hid_t sid = H5Screate_simple(rank, dims, NULL);

    /* Create the dataset. */
    hid_t did = H5Dcreate2(loc_id, dset_name, tid, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Write the dataset only if there is data to write */
    if (data && mpirank == 0)
    {
        status |= H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    }

    /* End access to the dataset and release resources used by it. */
    status |= H5Dclose(did);

    /* Terminate access to the data space. */
    status |= H5Sclose(sid);

    if (status != 0 || sid < 0 || did < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in make_single_value_dataset_numerical\n");
    }
}

/**
 * Computes chunk dimensions for compressed datasets.
 *
 * The target chunk payload is 8 MiB. Dimensions 1..rank-1 are either copied from
 * dataset_dims or set to 1 when limit_other_dims is enabled. The first dimension
 * is then reduced as needed to keep the total chunk payload within the target.
 *
 * @param tid The HDF5 datatype identifier.
 * @param rank The rank of the dataset.
 * @param dataset_dims The full dataset dimensions used as the initial chunk shape.
 * @param chunk_dims Output chunk dimensions.
 * @param limit_other_dims If nonzero, dimensions 1..rank-1 are forced to 1.
 */
static void
get_chunk_size(hid_t tid, 
               int rank, 
               const hsize_t *dataset_dims, 
               hsize_t *chunk_dims,
               int limit_other_dims)
{
    // target 8 MiB chunk size
    const hsize_t target_chunk_size = ((hsize_t) 8) << 20; /* 8 MiB */
    hsize_t chunk_size = H5Tget_size(tid);

    for(int i = 1; i < rank; i++)
    {
        chunk_dims[i] = limit_other_dims ? 1 : dataset_dims[i];

        if (chunk_size <= target_chunk_size && chunk_dims[i] > 0)
        {
            if (chunk_size > target_chunk_size / chunk_dims[i])
            {
                chunk_size = target_chunk_size + 1;
            }
            else
            {
                chunk_size *= chunk_dims[i];
            }
        }
    }

    hsize_t max_dim0 = 1;
    if (chunk_size > 0 && chunk_size <= target_chunk_size)
    {
        max_dim0 = target_chunk_size / chunk_size;
        if (max_dim0 < 1)
        {
            max_dim0 = 1;
        }
    }

    chunk_dims[0] = dataset_dims[0] < max_dim0 ? dataset_dims[0] : max_dim0;
    if (chunk_dims[0] < 1)
    {
        chunk_dims[0] = 1;
    }
}

/**
 * Creates a new HDF5 dataset. And returns the dataset identifier, 
 * so that the dataset can be written to.
 *
 * @param clawpatch_opts The clawpatch options.
 * @param loc_id The location identifier for the dataset.
 * @param tid The datatype for the dataset.
 * @param dset_name The name of the dataset.
 * @param rank The rank of the dataset.
 * @param dims The dimensions of the dataset.
 * @param patch_dims The dimensions for a single patch.
 * @return The dataset identifier.
 */
static hid_t 
make_dataset(const fclaw_clawpatch_options_t *clawpatch_opts,
             hid_t loc_id, 
             hid_t tid, 
             const char *dset_name, 
             int rank, 
             const hsize_t *dims, 
             const hsize_t *patch_dims)
{
    herr_t status = 0;

    hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);

    status |= H5Pset_fill_time(prop_id, H5D_FILL_TIME_NEVER);
    
    if(patch_dims != NULL && clawpatch_opts->hdf5_compression_level > 0)
    {
        hsize_t limited_chunk_dims[rank];
        get_chunk_size(tid,
                   rank,
                   patch_dims,
                   limited_chunk_dims,
                   0);
        status |= H5Pset_chunk(prop_id, rank, limited_chunk_dims);
        if(tid == H5T_NATIVE_INT || tid == H5T_NATIVE_UINT8)
        {
            status |= H5Pset_scaleoffset(prop_id, H5Z_SO_INT, 1);
        }
        else
        {
            status |= H5Pset_shuffle(prop_id);
        }
        status |= H5Pset_deflate(prop_id, clawpatch_opts->hdf5_compression_level);
    }

    
    /* Create the data space for the dataset. */
    hid_t sid = H5Screate_simple(rank, dims, NULL);

    /* Create the dataset. */
    hid_t did = H5Dcreate2(loc_id, dset_name, tid, sid, H5P_DEFAULT, prop_id, H5P_DEFAULT);

    status |= H5Sclose(sid);
    status |= H5Pclose(prop_id);

    if(status != 0 || sid < 0 || did < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in make_dataset\n");
    }

    return did;
}

static fclaw_clawpatch_vtk_vtable_entry_t* s_hdf5_vtk_entry = NULL;
static fclaw_vtk_cb_context_t s_hdf5_vtk_ctx;

static void
hdf5_vtk_entry_adapter(fclaw_global_t *glob,
                       fclaw_patch_t *patch,
                       int blockno,
                       int patchno,
                       char *buffer)
{
    FCLAW_ASSERT(s_hdf5_vtk_entry != NULL);
    s_hdf5_vtk_entry->callback(glob, patch, blockno, patchno, &s_hdf5_vtk_ctx, buffer);
}

static hid_t
hdf5_tid_from_vtk_type(fclaw_vtk_entry_type_t type, int fits32)
{
    switch (type)
    {
    case FCLAW_VTK_UINT8:
        return H5T_NATIVE_UINT8;
    case FCLAW_VTK_INT32:
        return H5T_NATIVE_INT32;
    case FCLAW_VTK_INT64:
        return H5T_NATIVE_INT64;
    case FCLAW_VTK_UINT64:
        return H5T_NATIVE_UINT64;
    case FCLAW_VTK_FLOAT32:
        return H5T_NATIVE_FLOAT;
    case FCLAW_VTK_FLOAT64:
        return H5T_NATIVE_DOUBLE;
    case FCLAW_VTK_INT32_OR_64:
        return fits32 ? H5T_NATIVE_INT32 : H5T_NATIVE_INT64;
    default:
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Unsupported vtk entry type\n");
        return H5T_NATIVE_INT32;
    }
}

static size_t
hdf5_entry_elements_in_patch(fclaw_global_t *glob,
                             fclaw_patch_t *patch,
                             int blockno,
                             int patchno,
                             const fclaw_clawpatch_vtk_vtable_entry_t *entry,
                             fclaw_vtk_cb_context_t *ctx)
{
    if (patch == NULL)
    {
        return 0;
    }
    return entry->elements_in_patch(glob, patch, blockno, patchno, ctx);
}

static hsize_t
write_vtable_entry_dataset(fclaw_global_t *glob,
                           hid_t loc_id,
                           const char *dset_name,
                           const fclaw_clawpatch_vtk_vtable_entry_t *entry,
                           int num_patches_to_buffer,
                           int fits32,
                           fclaw_hdf5_patch_data_t patch_cb)
{
    if (entry == NULL)
    {
        return 0;
    }

    const int rank = entry->number_of_components == 1 ? 1 : 2;
    hid_t tid = hdf5_tid_from_vtk_type(entry->type, fits32);

    fclaw_vtk_cb_context_t ctx;
    ctx.fits32 = fits32;
    ctx.offsets_include_zero = 1;

    unsigned long long local_total = 0;
    for (int local_patch_index = 0;
         local_patch_index < glob->domain->local_max_patches;
         ++local_patch_index)
    {
        int patchno, blockno;
        get_patch_blockno(glob, local_patch_index, &blockno, &patchno);
        fclaw_patch_t *patch = NULL;
        if (patchno < glob->domain->blocks[blockno].num_patches)
        {
            patch = &glob->domain->blocks[blockno].patches[patchno];
        }
        local_total += (unsigned long long)
            hdf5_entry_elements_in_patch(glob, patch, blockno, patchno, entry, &ctx);
    }

    unsigned long long global_start = 0;
    unsigned long long global_total = local_total;
#ifdef FCLAW_ENABLE_MPI
    unsigned long long prefix_total = 0;
    herr_t mpi_status = 0;
    mpi_status |= sc_MPI_Scan(&local_total, &prefix_total, 1,
                              sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM,
                              glob->mpicomm);
    mpi_status |= sc_MPI_Allreduce(&local_total, &global_total, 1,
                                   sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM,
                                   glob->mpicomm);
    if (mpi_status != 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in MPI scan/reduce\n");
    }
    global_start = prefix_total - local_total;
#endif

    hsize_t dataset_dims[2] = {0, 0};
    dataset_dims[0] = (hsize_t) global_total;
    if (rank == 2)
    {
        dataset_dims[1] = (hsize_t) entry->number_of_components;
    }

    hsize_t patch_dims[2] = {0, 0};
    patch_dims[0] = (hsize_t) entry->elements_per_patch;
    if (rank == 2)
    {
        patch_dims[1] = (hsize_t) entry->number_of_components;
    }

    hid_t did = make_dataset(fclaw_clawpatch_get_options(glob),
                             loc_id,
                             tid,
                             dset_name,
                             rank,
                             dataset_dims,
                             patch_dims);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    herr_t io_status = 0;

#ifdef FCLAW_ENABLE_MPI
    herr_t xfer_status = 0;
    xfer_status |= H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    if (xfer_status != 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in H5Pset_dxpl_mpio\n");
    }
#endif

    unsigned long long local_prefix = 0;
    for (int local_patch_index = 0;
         local_patch_index < glob->domain->local_max_patches;
         local_patch_index += num_patches_to_buffer)
    {
        int num_able_to_buffer = SC_MIN(num_patches_to_buffer,
                                        glob->domain->local_max_patches - local_patch_index);

        hid_t filespace = H5Dget_space(did);
        hsize_t slab_dims[2] = {0, 0};
        hsize_t slab_start[2] = {0, 0};
        hsize_t buffer_offsets[num_able_to_buffer];
        hsize_t total_buffer_size = 0;
        hsize_t batch_total = 0;

        for (int j = 0; j < num_able_to_buffer; j++)
        {
            int patchno, blockno;
            get_patch_blockno(glob, local_patch_index + j, &blockno, &patchno);

            fclaw_patch_t *patch = NULL;
            if (patchno < glob->domain->blocks[blockno].num_patches)
            {
                patch = &glob->domain->blocks[blockno].patches[patchno];
            }

            size_t patch_count = hdf5_entry_elements_in_patch(glob, patch, blockno, patchno, entry, &ctx);
            buffer_offsets[j] = total_buffer_size;
            hsize_t buffer_size = (hsize_t) patch_count * (hsize_t) entry->number_of_components * (hsize_t) H5Tget_size(tid);
            total_buffer_size += buffer_size;
            batch_total += (hsize_t) patch_count;

            if (j == 0)
            {
                slab_start[0] = (hsize_t) global_start + (hsize_t) local_prefix;
                slab_dims[0] = (hsize_t) patch_count;
                if (rank == 2)
                {
                    slab_start[1] = 0;
                    slab_dims[1] = (hsize_t) entry->number_of_components;
                }
            }
            else
            {
                slab_dims[0] += (hsize_t) patch_count;
            }
        }

        local_prefix += batch_total;

        io_status = 0;
        io_status |= H5Sselect_hyperslab(filespace, H5S_SELECT_SET, slab_start, NULL, slab_dims, NULL);
        hid_t memspace = H5Screate_simple(rank, slab_dims, NULL);
        FCLAW_ASSERT(H5Sget_select_npoints(filespace) == H5Sget_select_npoints(memspace));

        char *buffer = FCLAW_ALLOC(char, total_buffer_size);
        for (int j = 0; j < num_able_to_buffer; j++)
        {
            int patchno, blockno;
            get_patch_blockno(glob, local_patch_index + j, &blockno, &patchno);
            if (buffer_offsets[j] < total_buffer_size)
            {
                fclaw_patch_t *patch = NULL;
                if (patchno < glob->domain->blocks[blockno].num_patches)
                {
                    patch = &glob->domain->blocks[blockno].patches[patchno];
                }
                if (patch_cb != NULL)
                {
                    patch_cb(glob, patch, blockno, patchno, buffer + buffer_offsets[j]);
                }
                else
                {
                    s_hdf5_vtk_entry = (fclaw_clawpatch_vtk_vtable_entry_t *) entry;
                    s_hdf5_vtk_ctx = ctx;
                    hdf5_vtk_entry_adapter(glob, patch, blockno, patchno, buffer + buffer_offsets[j]);
                }
            }
        }

        io_status |= H5Dwrite(did, tid, memspace, filespace, plist_id, buffer);
        io_status |= H5Sclose(memspace);
        io_status |= H5Sclose(filespace);
        FCLAW_FREE(buffer);
        if (io_status != 0)
        {
            fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in write_vtable_entry_dataset\n");
        }
    }

    io_status = 0;
    io_status |= H5Pclose(plist_id);
    io_status |= H5Dclose(did);
    if (io_status != 0 || did < 0 || plist_id < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in write_vtable_entry_dataset\n");
    }

    return (hsize_t) global_total;
}

//this is copied form hdf5lt needed to change to NULLPAD for VTK
void
set_attribute_string(hid_t loc_id, const char *obj_name, const char *attr_name, const char *attr_data)
{
    herr_t status = 0;

    /* Open the object */
    hid_t obj_id = H5Oopen(loc_id, obj_name, H5P_DEFAULT);

    hid_t attr_type = H5Tcopy(H5T_C_S1);

    hsize_t attr_size = strlen(attr_data);

    status |= H5Tset_size(attr_type, (size_t)attr_size);
    status |= H5Tset_strpad(attr_type, H5T_STR_NULLPAD);

    hid_t attr_space_id = H5Screate(H5S_SCALAR);

    /* Delete the attribute if it already exists */
    hid_t has_attr = H5Aexists(obj_id, attr_name);
    if (has_attr > 0)
    {
        status |= H5Adelete(obj_id, attr_name);
    }

    /* Create and write the attribute */

    hid_t attr_id = H5Acreate2(obj_id, attr_name, attr_type, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);

    status |= H5Awrite(attr_id, attr_type, attr_data);
    status |= H5Aclose(attr_id);
    status |= H5Sclose(attr_space_id);
    status |= H5Tclose(attr_type);
    status |= H5Oclose(obj_id);

    if(status != 0 || obj_id < 0 || attr_type < 0 || attr_space_id < 0 || attr_id < 0)
    {
        printf("fclaw_clawpatch_output_hdf5.c Error in set_attribute_string\n");
    }
}

/**
 * @brief Write vtkhdf file
 * 
 * @param glob the global context
 * @param filename the name of the file
 * @param coordinate_cb the callback function to get the coordinate data
 * @param value_cb the callback function to get the value data
 */
static hsize_t
hdf5_entry_global_count(fclaw_global_t *glob,
                        const fclaw_clawpatch_vtk_vtable_entry_t *entry,
                        fclaw_vtk_cb_context_t *ctx)
{
    if (entry == NULL)
    {
        return 0;
    }

    unsigned long long local_total = 0;
    for (int local_patch_index = 0;
         local_patch_index < glob->domain->local_max_patches;
         ++local_patch_index)
    {
        int patchno, blockno;
        get_patch_blockno(glob, local_patch_index, &blockno, &patchno);
        fclaw_patch_t *patch = NULL;
        if (patchno < glob->domain->blocks[blockno].num_patches)
        {
            patch = &glob->domain->blocks[blockno].patches[patchno];
        }
        local_total += (unsigned long long)
            hdf5_entry_elements_in_patch(glob, patch, blockno, patchno, entry, ctx);
    }

    unsigned long long global_total = local_total;
#ifdef FCLAW_ENABLE_MPI
    unsigned long long prefix_total = 0;
    herr_t mpi_status = 0;
    mpi_status |= sc_MPI_Scan(&local_total, &prefix_total, 1,
                              sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM,
                              glob->mpicomm);
    mpi_status |= sc_MPI_Allreduce(&local_total, &global_total, 1,
                                   sc_MPI_UNSIGNED_LONG_LONG, sc_MPI_SUM,
                                   glob->mpicomm);
    if (mpi_status != 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in MPI scan/reduce\n");
    }
#endif

    return (hsize_t) global_total;
}

static void
fclaw_hdf_write_file (fclaw_global_t * glob, 
                      const char* filename,
                      fclaw_hdf5_patch_data_t coordinate_cb,
                      fclaw_hdf5_patch_data_t value_cb)
{
    const fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);
    fclaw_clawpatch_vtk_vtable_t *vtk_vtable = fclaw_clawpatch_vtk_vtable(glob);

    int num_patches_to_buffer = clawpatch_opt->hdf5_patch_threshold;
    if(clawpatch_opt->hdf5_patch_threshold == 0)
    {
        num_patches_to_buffer = glob->domain->local_max_patches;
    }

    char vtkhdf[8] = "/VTKHDF";
    char celldata[18] = "/VTKHDF/CellData";
    
    fclaw_remove(filename);
    herr_t status = 0;
    // Set up file access property list with parallel I/O access
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef FCLAW_ENABLE_MPI
    status |= H5Pset_fapl_mpio(fapl_id, glob->mpicomm, MPI_INFO_NULL);
    status |= H5Pset_coll_metadata_write(fapl_id, 1);
    status |= H5Pset_all_coll_metadata_ops(fapl_id, 1);
#endif
    status |= H5Pset_libver_bounds(fapl_id, H5F_LIBVER_V110, H5F_LIBVER_LATEST);

    hid_t fcpl_id = H5Pcreate(H5P_FILE_CREATE);

    // Create a new file collectively and release property list identifier.
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, fcpl_id, fapl_id);

    status |= H5Pclose(fapl_id);
    status |= H5Pclose(fcpl_id);

    hid_t vtkhdf_gid = H5Gcreate2(file_id, vtkhdf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    int vtk_version[2] = {1, 0};
    set_attribute_numerical(file_id, vtkhdf, "Version", 2, H5T_NATIVE_INT, vtk_version);
    set_attribute_string(file_id, vtkhdf, "Type", "UnstructuredGrid");
    
    // The idea here is to write a paritioned vtu, where each partition is a patch
    // this will make processing easier in matlab

    hsize_t dims[3] = {0,0,0};

    // write single value datasets for vtk

    fclaw_vtk_cb_context_t size_ctx;
    size_ctx.fits32 = 0;
    size_ctx.offsets_include_zero = 1;

    hsize_t number_of_points = hdf5_entry_global_count(glob, &vtk_vtable->position_entry, &size_ctx);
    hsize_t number_of_connectivity_ids = hdf5_entry_global_count(glob, &vtk_vtable->connectivity_entry, &size_ctx);
    hsize_t number_of_cells = hdf5_entry_global_count(glob, &vtk_vtable->types_entry, &size_ctx);

    int fits32 = number_of_points <= INT32_MAX
        && number_of_connectivity_ids <= INT32_MAX
        && number_of_cells <= INT32_MAX;

    s_hdf5_vtk_ctx.fits32 = fits32;
    hsize_t written_number_of_points = write_vtable_entry_dataset(glob,
                                                                  vtkhdf_gid,
                                                                  "Points",
                                                                  &vtk_vtable->position_entry,
                                                                  num_patches_to_buffer,
                                                                  fits32,
                                                                  coordinate_cb);
    hsize_t written_number_of_connectivity_ids = write_vtable_entry_dataset(glob,
                                                                            vtkhdf_gid,
                                                                            "Connectivity",
                                                                            &vtk_vtable->connectivity_entry,
                                                                            num_patches_to_buffer,
                                                                            fits32,
                                                                            NULL);
    hsize_t written_number_of_cells = write_vtable_entry_dataset(glob,
                                                                 vtkhdf_gid,
                                                                 "Offsets",
                                                                 &vtk_vtable->offsets_entry,
                                                                 num_patches_to_buffer,
                                                                 fits32,
                                                                 NULL);
    write_vtable_entry_dataset(glob,
                               vtkhdf_gid,
                               "Types",
                               &vtk_vtable->types_entry,
                               num_patches_to_buffer,
                               fits32,
                               NULL);

    FCLAW_ASSERT(written_number_of_points == number_of_points);
    FCLAW_ASSERT(written_number_of_connectivity_ids == number_of_connectivity_ids);
    FCLAW_ASSERT(written_number_of_cells == number_of_cells + 1);

    dims[0] = 1;
    long number_of_cells_long = (long) number_of_cells;
    make_single_value_dataset_numerical(glob->mpirank, vtkhdf_gid, "NumberOfCells", 1, dims, H5T_NATIVE_LONG, &number_of_cells_long);

    long number_of_points_long = (long) number_of_points;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, vtkhdf_gid, "NumberOfPoints", 1, dims, H5T_NATIVE_LONG, &number_of_points_long);

    long number_of_connectivity_ids_long = (long) number_of_connectivity_ids;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, vtkhdf_gid, "NumberOfConnectivityIds", 1, dims, H5T_NATIVE_LONG, &number_of_connectivity_ids_long);

    hid_t fielddata_gid = H5Gcreate2(file_id,
                                     "/VTKHDF/FieldData",
                                     H5P_DEFAULT,
                                     H5P_DEFAULT,
                                     H5P_DEFAULT);
    sc_link_t* curr_entry = vtk_vtable->field_entries->first;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        write_vtable_entry_dataset(glob,
                                   fielddata_gid,
                                   entry->name,
                                   entry,
                                   num_patches_to_buffer,
                                   fits32,
                                   NULL);
        curr_entry = curr_entry->next;
    }
    status |= H5Gclose(fielddata_gid);

    /* avoid resource leaks by closing */
    status |= H5Gclose(vtkhdf_gid);


    hid_t celldata_gid = H5Gcreate2(file_id, celldata, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    curr_entry = vtk_vtable->celldata_entries->first;
    while (curr_entry != NULL)
    {
        fclaw_clawpatch_vtk_vtable_entry_t* entry =
            (fclaw_clawpatch_vtk_vtable_entry_t*) curr_entry->data;
        fclaw_hdf5_patch_data_t patch_cb = NULL;
        if (strcmp(entry->name, "meqn") == 0)
        {
            patch_cb = value_cb;
        }
        write_vtable_entry_dataset(glob,
                                   celldata_gid,
                                   entry->name,
                                   entry,
                                   num_patches_to_buffer,
                                   fits32,
                                   patch_cb);
        curr_entry = curr_entry->next;
    }

    status |= H5Gclose(celldata_gid);

    status |= H5Fclose(file_id);
    
    FCLAW_ASSERT(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL) == 0);

    if(status != 0 || file_id < 0 || vtkhdf_gid < 0 || celldata_gid < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in fclaw_hdf_write_file\n");
    }
}

/*  ---------------------------------------------------------------------------
    Public interface
    --------------------------------------------------------------------------- */

void fclaw_clawpatch_output_hdf5_to_file (struct fclaw_global* glob, 
                                         const char* filename,
                                         fclaw_hdf5_patch_data_t coordinate_cb,
                                         fclaw_hdf5_patch_data_t value_cb)
{
    fclaw_hdf_write_file (glob, 
                          filename, 
                          coordinate_cb,
                          value_cb);
}

void fclaw_clawpatch_output_hdf5 (fclaw_global_t * glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d.vtkhdf", fclaw_opt->prefix, iframe);

    fclaw_clawpatch_output_hdf5_to_file(glob,basename, NULL, NULL);
}



