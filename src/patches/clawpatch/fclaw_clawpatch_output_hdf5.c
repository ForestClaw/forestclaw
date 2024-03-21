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
#include <fclaw_map.h>
#include <hdf5.h>

/*----------------------------------------------------------------------
    Value Callbacks
  ----------------------------------------------------------------------*/
static int get_num_cells(fclaw_global_t * glob,
                         fclaw_patch_t* patch)
{
    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;

    int num_cells;
    if(fclaw_clawpatch_dim(patch) == 2)
    {
        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,&xlower,&ylower,&dx,&dy);
        num_cells = mx * my;
    }
    else
    {
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz,&mbc,&xlower,&ylower,&zlower,&dx,&dy,&dz);
        num_cells = mx * my * mz;
    }
    return num_cells;
}

static void 
blockno_cb (fclaw_global_t * glob,
            fclaw_patch_t * patch,
            int blockno,
            int patchno,
            char * buffer)
{
    int *blocknos = (int *) buffer;

    int num_cells = get_num_cells(glob,patch);

    for(int i = 0; i < num_cells; i++)
    {
        blocknos[i] = blockno;
    }
}

static void 
patchno_cb (fclaw_global_t * glob,
            fclaw_patch_t * patch,
            int blockno,
            int patchno,
            char * buffer)
{
    int *patchnos = (int *) buffer;

    int num_cells = get_num_cells(glob,patch);

    for(int i = 0; i < num_cells; i++)
    {
        patchnos[i] = patchno;
    }
}

static void 
mpirank_cb (fclaw_global_t * glob,
            fclaw_patch_t * patch,
            int blockno,
            int patchno,
            char * buffer)
{
    int *mpiranks = (int *) buffer;

    int num_cells = get_num_cells(glob,patch);

    for(int i = 0; i < num_cells; i++)
    {
        mpiranks[i] = glob->mpirank;
    }
}

static void 
types_cb (fclaw_global_t * glob,
          fclaw_patch_t * patch,
          int blockno,
          int patchno,
          char * buffer)
{
    uint8_t * types = (uint8_t *) buffer;

    int num_cells = get_num_cells(glob,patch);

    uint8_t type;
    if(fclaw_clawpatch_dim(patch) == 2)
    {
        type = 9;
    }
    else
    {
        type = 12;
    }

    for(int i = 0; i < num_cells; i++)
    {
        types[i] = type;
    }
}

static void
get_coordinates (fclaw_global_t * glob,
                 fclaw_patch_t * patch,
                 int blockno,
                 int patchno,
                 char * buffer)
{
    double * points = (double *) buffer;
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    const fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int mx,my,mz,mbc;
    double dx,dy,dz,xlower,ylower,zlower;

    if(clawpatch_opt->patch_dim == 2)
    {
        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);
        mz = 0;
        dz = 0;
        zlower = 0;
    }
    else 
    {
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc,
                                    &xlower,&ylower,&zlower, &dx,&dy, &dz);
    }

    fclaw_map_context_t* cont = fclaw_map_get(glob);
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
                if (clawpatch_opt->patch_dim == 2 && fclaw_opt->manifold)
                {
                    FCLAW_MAP_2D_C2M(&cont,&blockno,&x,&y,&xpp,&ypp,&zpp);
                    *points++ = xpp;
                    *points++ = ypp;
                    *points++ = zpp;
                }
                else if (clawpatch_opt->patch_dim == 3 && fclaw_opt->manifold)
                {
                    FCLAW_MAP_3D_C2M(&cont,&blockno,&x,&y,&z,&xpp,&ypp,&zpp);
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

#define GET_OFFSETS(TYPE) \
static void \
get_offsets_##TYPE(fclaw_global_t * glob, \
            fclaw_patch_t * patch, \
            int blockno, \
            int patchno, \
            char * buffer) \
{ \
    TYPE * offsets = (TYPE *) buffer; \
\
    int num_cells_before; \
    int mx,my,mz,mbc; \
    double dx,dy,dz,xlower,ylower,zlower; \
\
    int num_patches_before = glob->domain->global_num_patches_before + glob->domain->blocks[blockno].num_patches_before + patchno; \
    int num_points_per_cell; \
    if(fclaw_clawpatch_dim(patch) == 2) \
    { \
        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,&xlower,&ylower,&dx,&dy); \
        mz = 1; \
        num_points_per_cell = 4; \
        num_cells_before = num_patches_before * mx * my; \
    } \
    else \
    { \
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz,&mbc,&xlower,&ylower,&zlower,&dx,&dy,&dz); \
        num_points_per_cell = 8; \
        num_cells_before = num_patches_before * mx * my * mz; \
    } \
\
    int curr_offset = num_cells_before * num_points_per_cell; \
\
    for(int k = 0; k < mz; k++) \
    { \
        for(int j = 0; j < my; j++) \
        { \
            for(int i = 0; i < mx; i++) \
            { \
                *offsets++ = curr_offset; \
                curr_offset += num_points_per_cell; \
            } \
        } \
    } \
\
    if(glob->mpirank == glob->mpisize - 1 \
       && blockno == glob->domain->num_blocks - 1 \
       && patchno == glob->domain->blocks[blockno].num_patches - 1) \
    { \
        /* last patch on last rank has extra value */ \
       *offsets = curr_offset; \
    } \
}

GET_OFFSETS(int32_t)
GET_OFFSETS(int64_t)

#define POINT_INDEX_CB(TYPE) \
static void \
get_point_indexes_##TYPE (fclaw_global_t * glob, \
                          fclaw_patch_t * patch, \
                          int blockno, \
                          int patchno, \
                          char * buffer) \
{ \
    TYPE * points = (TYPE *) buffer; \
    const fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob); \
\
    int mx,my,mz,mbc; \
    double dx,dy,dz,xlower,ylower,zlower; \
\
    TYPE num_points_before = glob->domain->global_num_patches_before + glob->domain->blocks[blockno].num_patches_before + patchno; \
    if(clawpatch_opt->patch_dim == 2) \
    { \
        fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc, \
                                    &xlower,&ylower,&dx,&dy); \
        mz = 0; \
        dz = 0; \
        zlower = 0; \
        num_points_before *= (mx + 1) * (my + 1); \
    } \
    else \
    { \
        fclaw_clawpatch_3d_grid_data(glob,patch,&mx,&my,&mz, &mbc, \
                                    &xlower,&ylower,&zlower, &dx,&dy, &dz); \
        num_points_before *= (mx + 1) * (my + 1) * (mz + 1); \
    } \
\
    /* Enumerate point coordinates in the patch */ \
    int i, j, k; \
    for (k = 0; k <= mz; ++k) \
    { \
        for (j = 0; j <= my; ++j) \
        { \
            for (i = 0; i <= mx; ++i) \
            { \
                *points = num_points_before; \
                points++; \
                num_points_before++; \
            } \
        } \
    } \
}

POINT_INDEX_CB(int32_t)
POINT_INDEX_CB(int64_t)

static void
write_2d_patch_q (fclaw_global_t * glob,
                  fclaw_patch_t * patch,
                  int blockno,
                  int patchno,
                  float * q_out)
{
    int mx,my,mbc;
    double dx,dy,xlower,ylower;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my, &mbc,
                                &xlower,&ylower, &dx,&dy);

    int meqn;
    double* q;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    for (int j = 0; j < my; ++j)
    {
        for (int i = 0; i < mx; ++i)
        {
            for (int m = 0; m < meqn; m++)
            {
                size_t out_idx = j*mx*meqn + i*meqn + m;
                size_t in_idx = m*(mx+2*mbc)*(my+2*mbc) + (j+mbc)*(mx+2*mbc) + (i+mbc);
                q_out[out_idx] = q[in_idx];
            }
        }
    }
}
static void
write_3d_patch_q (fclaw_global_t * glob,
                  fclaw_patch_t * patch,
                  int blockno,
                  int patchno,
                  float * q_out)
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
get_data (fclaw_global_t * glob,
          fclaw_patch_t * patch,
          int blockno,
          int patchno,
          char * buffer)
{
    float * q_out = (float *) buffer;
    if (fclaw_clawpatch_dim(patch) == 2)
    {
        write_2d_patch_q(glob, patch, blockno, patchno, q_out);
    }
    else
    {
        write_3d_patch_q(glob, patch, blockno, patchno, q_out);
    }
}

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
        blockno++;
    }
    *patchno = i - glob->domain->blocks[*blockno].num_patches_before;
}

/**
 * Check if the given block number and patch number are valid.
 *
 * @param glob The global context.
 * @param blockno The block number to check.
 * @param patchno The patch number to check.
 * @return 1 if the block number and patch number are valid, 0 otherwise.
 */
static int 
patch_blockno_valid(fclaw_global_t *glob, int blockno, int patchno)
{
    return blockno < glob->domain->num_blocks && patchno < glob->domain->blocks[blockno].num_patches;
}

/*----------------------------------------------------------------------
    Virtualized funcitons
  ----------------------------------------------------------------------*/

/**
 * @brief Virtualized functions needed for making a patch dataset
 */
typedef struct make_dataset_vtable
{
    /**
     * @brief Get the overall dimensions of a dataset
     *
     * @param[in] glob The global context.
     * @param[in] rank The rank of the patch.
     * @param[in] patch_dims The dimensions of the patch.
     * @param[out] dataset_dims The array to store the dataset dimensions.
     */
    void (*get_dataset_dims)(fclaw_global_t *glob, 
                             int rank, 
                             const hsize_t *patch_dims, 
                             hsize_t *dataset_dims);
    /**
     * @brief Get the slab dimensions for a given patch in a block.
     *
     * @param[in] glob The global context.
     * @param[in] patchno The patch number.
     * @param[in] blockno The block number.
     * @param[in] rank The rank of the data.
     * @param[in] patch_dims The dimensions of the patch.
     * @param[out] slab_start The starting indices of the slab.
     * @param[out] slab_dims The dimensions of the slab.
     */
    void (*get_slab_dims)(fclaw_global_t *glob, 
                          int patchno, 
                          int blockno, int rank, 
                          const hsize_t *patch_dims, 
                          hsize_t *slab_start, 
                          hsize_t *slab_dims);
} make_dataset_vtable_t;

/* ------ default functions -----*/

/**
 * Calculates the slab dimensions for a given patch in a block.
 * This just returns the patch dimensions, and the starting index.
 * If the patch is not valid, the slab dimensions are set to 0.
 *
 * @param[in] glob The global context.
 * @param[in] patchno The patch number.
 * @param[in] blockno The block number.
 * @param[in] rank The rank of the data.
 * @param[in] patch_dims The dimensions of the patch.
 * @param[out] slab_start The starting indices of the slab.
 * @param[out] slab_dims The dimensions of the slab.
 */
static void
get_slab_dims(fclaw_global_t *glob, 
              int patchno, int blockno, 
              int rank, const hsize_t *patch_dims, hsize_t *slab_start, hsize_t *slab_dims)
{
    if(patch_blockno_valid(glob, blockno, patchno))
    {
        int i = glob->domain->blocks[blockno].num_patches_before + patchno;
        slab_start[0] = (glob->domain->global_num_patches_before + i) * patch_dims[0];
        slab_dims[0] = patch_dims[0];
        for(int i = 1; i < rank; i++)
        {
            slab_start[i] = 0;
            slab_dims[i] = patch_dims[i];
        }
    }
    else
    {
        for(int i = 0; i < rank; i++)
        {
            slab_start[i] = 0;
            slab_dims[i] = 0;
        }
    }
}

/**
 * @brief Calculate the overall dataset dimensions.
 * This is just the first patch dimensions multiplied by the number of patches.
 *
 * @param[in] glob The global context.
 * @param[in] rank The rank of the patch.
 * @param[in] patch_dims The dimensions of the patch.
 * @param[out] dataset_dims The array to store the dataset dimensions.
 */
static 
void get_dataset_dims(fclaw_global_t *glob, 
                      int rank, 
                      const hsize_t *patch_dims, 
                      hsize_t *dataset_dims)
{
    for(int i = 0; i < rank; i++)
    {
        dataset_dims[i] = patch_dims[i];
    }
    dataset_dims[0] *= glob->domain->global_num_patches;
}

/** default vtable */
static make_dataset_vtable_t default_vtable = 
{
    get_dataset_dims,
    get_slab_dims
};
/**
 * @brief Calculate the overall dataset dimensions for the offset dataset
 * This is the same as default plus an extra value for the last patch on the last rank.
 *
 * @param[in] glob The global context.
 * @param[in] rank The rank of the patch.
 * @param[in] patch_dims The dimensions of the patch.
 * @param[out] dataset_dims The array to store the dataset dimensions.
 */
static void
get_dataset_dims_offset(fclaw_global_t *glob, int rank, const hsize_t *patch_dims, hsize_t *dataset_dims)
{
    get_dataset_dims(glob, rank, patch_dims, dataset_dims);
    dataset_dims[0] += 1;
}
/**
 * Same as default but on last patch on last rank, the slab dimensions are increased by 1.
 *
 * @param[in] glob The global context.
 * @param[in] patchno The patch number.
 * @param[in] blockno The block number.
 * @param[in] rank The rank of the data.
 * @param[in] patch_dims The dimensions of the patch.
 * @param[out] slab_start The starting indices of the slab.
 * @param[out] slab_dims The dimensions of the slab.
 */
static void
get_slab_dims_offset(fclaw_global_t *glob, 
                     int patchno, int blockno, 
                     int rank, const hsize_t *patch_dims, hsize_t *slab_start, hsize_t *slab_dims)
{
    get_slab_dims(glob, patchno, blockno, rank, patch_dims, slab_start, slab_dims);
    if(glob->mpirank == glob->mpisize - 1 
       && blockno == glob->domain->num_blocks - 1 
       && patchno == glob->domain->blocks[blockno].num_patches - 1)
    {
        // last patch on last rank has extra value
        slab_dims[0] += 1;
    }
}
/** vtable for offset datset */
static make_dataset_vtable_t offset_vtable = 
{
    get_dataset_dims_offset,
    get_slab_dims_offset
};

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
 * Limits the chunk size for HDF5.
 * The initial chunk size is a contiguous block of patch data.
 * The max limit of a chunk in HDF5 is 4GB. So there can be cases where 
 * the chunk size has to be limited to fit within the 4GB limit.
 *
 * @param tid The HDF5 datatype identifier.
 * @param subdims_rank The rank of the subdimensions.
 * @param subdims The subdimensions array. 
                  The first dimenion of chunk_dims is a flattend array of patch data.
                  For example if there is a 2D patch data that is 8x8 cells, then
                  chunk_dims[0] will be 64 and the subdims array will be {8,8}.
 * @param rank The rank of the data
 * @param chunk_dims The chunk dimensions array.
 * @param limited_chunk_dims The limited chunk dimensions array.
 * @param limit_other_dims Flag indicating whether to limit other dimensions for better compression.
                           The data across the other dimensions will be dissimilar. So the 
                           other chunk dimensions can be set to 1 in order to compress better, 
                           since the chunks will have similar data.
 */
static void
limit_chunk_size(hid_t tid, 
                 int subdims_rank,
                 const hsize_t *subdims,
                 int rank, 
                 const hsize_t *chunk_dims, 
                 hsize_t *limited_chunk_dims,
                 int limit_other_dims)
{
    hsize_t max_chunk_size = 1 << 28; // 4*2^30 bytes, 4GB
    hsize_t type_size = H5Tget_size(tid);
    hsize_t chunk_size = type_size;
    if(limit_other_dims)
    {
        //only first dimension is similar data,
        //so chunk other dimensions for better compression
        for(int i = 1; i < rank; i++)
        {
            limited_chunk_dims[i] = 1;
        }
    }
    else
    {
        for(int i = rank-1; i >= 1; i--)
        {
            chunk_size *= chunk_dims[i];
            if(chunk_size > max_chunk_size)
            {
                limited_chunk_dims[i] = 1;
            }
            else
            {
                limited_chunk_dims[i] = chunk_dims[i];
            }
        }
    }
    // the subdims in forestclaw output are column-major
    // check that we aren't exceeding 4gb
    if(subdims_rank == 0)
    {
        limited_chunk_dims[0] = 1;
    }
    else
    {
        limited_chunk_dims[0] = 1;
        for(int i = 0; i < subdims_rank; i++)
        {
            chunk_size *= subdims[i];
            if(chunk_size > max_chunk_size)
            {
                break;
            }
            else
            {
                limited_chunk_dims[0] *= subdims[i];
            }
        }
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
 * @param subdims_rank The rank of the sub-dimensions.
 * @param subdims The sub-dimensions of the dataset. See limit_chunk_size for more information.
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
             int subdims_rank,
             const hsize_t *subdims,
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
        int limit_other_dims = tid != H5T_NATIVE_INT32 && tid != H5T_NATIVE_INT64;
        limit_chunk_size(tid, 
                         subdims_rank, subdims,
                         rank, patch_dims, limited_chunk_dims,
                         limit_other_dims);
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

    if(status != 0 || sid < 0 || did < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in make_dataset\n");
    }

    return did;
}

/**
 * @brief Makes a patch dataset.
 * 
 * @param glob the global context
 * @param loc_id the location identifier for the dataset
 * @param dset_name the name of the dataset
 * @param subdims_rank the rank of the sub-dimensions
 * @param subdims the sub-dimensions of the dataset. see limit_chunk_size for more information
 * @param rank the rank of the dataset
 * @param patch_dims the dimensions for a single patch
 * @param num_patches_to_buffer the number of patches to buffer
 * @param tid the datatype for the dataset 
 * @param patch_cb the callback function to get the patch data
 * @param vt the virtualized functions for creating the dataset
 */
static void
make_dataset_numerical(fclaw_global_t *glob,
                       hid_t loc_id, 
                       const char *dset_name, 
                       int subdims_rank, 
                       const hsize_t *subdims,
                       int rank, 
                       const hsize_t *patch_dims, 
                       int num_patches_to_buffer,
                       hid_t tid,
                       fclaw_hdf5_patch_data_t patch_cb,
                       make_dataset_vtable_t *vt)
{
    const fclaw_clawpatch_options_t *clawpatch_opts = fclaw_clawpatch_get_options(glob);

    herr_t status = 0;

    hsize_t dataset_dims[rank];
    vt->get_dataset_dims(glob, rank, patch_dims, dataset_dims);

    /* Create the dataset. */
    hid_t did = make_dataset(clawpatch_opts, 
                             loc_id, 
                             tid, 
                             dset_name, 
                             subdims_rank, subdims,
                             rank, dataset_dims, patch_dims);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

#ifdef FCLAW_ENABLE_MPI
    status |= H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    //status |= H5Pset_dxpl_mpio_collective_opt(plist_id, H5FD_MPIO_INDIVIDUAL_IO);
#endif

    for(int local_patch_index = 0; 
        local_patch_index < glob->domain->local_max_patches; 
        local_patch_index += num_patches_to_buffer)
    {
        int num_able_to_buffer 
            = SC_MIN(num_patches_to_buffer, glob->domain->local_max_patches - local_patch_index);

        hid_t filespace = H5Dget_space(did);

        hsize_t slab_dims[rank];
        hsize_t slab_start[rank];
        hsize_t buffer_offsets[num_able_to_buffer];
        hsize_t total_buffer_size = 0;
        for(int j=0; j < num_able_to_buffer; j++)
        {
            int patchno, blockno;
            get_patch_blockno(glob, local_patch_index+j, &blockno, &patchno);

            hsize_t patch_slab_dims[rank];
            hsize_t patch_slab_start[rank];
            vt->get_slab_dims(glob, patchno, blockno, rank, patch_dims, patch_slab_start, patch_slab_dims);
            buffer_offsets[j] = total_buffer_size;
            hsize_t buffer_size = patch_slab_dims[0] * H5Tget_size(tid);
            for(int k = 1; k < rank; k++)
            {
                buffer_size *= patch_slab_dims[k];
            }
            total_buffer_size += buffer_size;
            if(j == 0)
            {
                for(int k = 0; k < rank; k++)
                {
                    slab_dims[k] = patch_slab_dims[k];
                    slab_start[k] = patch_slab_start[k];
                }
            }
            else
            {
                slab_dims[0] += patch_slab_dims[0];
            }
        }

        status |= H5Sselect_hyperslab(filespace, H5S_SELECT_SET, slab_start, NULL, slab_dims, NULL);
        //create memspace
        hid_t memspace = H5Screate_simple(rank, slab_dims, NULL);
        FCLAW_ASSERT(H5Sget_select_npoints(filespace) == H5Sget_select_npoints(memspace));

        char *buffer = FCLAW_ALLOC(char, total_buffer_size);
        for(int j=0; j < num_able_to_buffer; j++)
        {
            int patchno, blockno;
            get_patch_blockno(glob, local_patch_index+j, &blockno, &patchno);
            if(buffer_offsets[j] < total_buffer_size)
            {
                fclaw_patch_t *patch = NULL;
                if(patchno < glob->domain->blocks[blockno].num_patches)
                {
                    patch = &glob->domain->blocks[blockno].patches[patchno];
                }
                patch_cb(glob, patch, blockno, patchno, buffer + buffer_offsets[j]);
            }
        }

        status |= H5Dwrite(did, tid, memspace, filespace, plist_id, buffer);

        status |= H5Sclose(memspace);
        status |= H5Sclose(filespace);

        FCLAW_FREE(buffer);
    }

    status |= H5Pclose(plist_id);
    status |= H5Dclose(did);

    if(status != 0 || did < 0 || plist_id < 0)
    {
        fclaw_abortf("fclaw_clawpatch_output_hdf5.c Error in make_dataset_numerical\n");
    }
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
static void
fclaw_hdf_write_file (fclaw_global_t * glob, 
                      const char* filename,
                      fclaw_hdf5_patch_data_t coordinate_cb,
                      fclaw_hdf5_patch_data_t value_cb)
{
    const fclaw_clawpatch_options_t* clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int num_patches_to_buffer = clawpatch_opt->hdf5_patch_threshold;
    if(clawpatch_opt->hdf5_patch_threshold == 0)
    {
        num_patches_to_buffer = glob->domain->local_max_patches;
    }

    //get mx, my, mz, meqn from clawpatch options
    int patch_dim = clawpatch_opt->patch_dim;
    int mx   = clawpatch_opt->mx;
    int my   = clawpatch_opt->my;
    int mz   = clawpatch_opt->mz;

    int global_num_patches = glob->domain->global_num_patches;

    int num_cells_per_patch;
    int num_points_per_patch;
    int num_points_per_cell;
    if(clawpatch_opt->patch_dim == 2)
    {
        num_cells_per_patch = mx * my;
        num_points_per_patch = (mx + 1) * (my + 1);
        num_points_per_cell = 4;
    }
    else 
    {
        num_cells_per_patch = mx * my * mz;
        num_points_per_patch = (mx + 1) * (my + 1) * (mz + 1);
        num_points_per_cell = 8;
    }
    hsize_t  num_cells_per_patch_subdims[3] = {mx, my, mz};
    hsize_t  num_points_per_patch_subdims[3] = {mx+1, my+1, mz+1};

    char vtkhdf[8] = "/VTKHDF";
    char celldata[18] = "/VTKHDF/CellData";
    
    herr_t status = 0;
    // Set up file access property list with parallel I/O access
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef FCLAW_ENABLE_MPI
    status |= H5Pset_fapl_mpio(fapl_id, glob->mpicomm, MPI_INFO_NULL);
    status |= H5Pset_coll_metadata_write(fapl_id, 1);
    status |= H5Pset_all_coll_metadata_ops(fapl_id, 1);
#endif

    hid_t fcpl_id = H5Pcreate(H5P_FILE_CREATE);
    status |= H5Pset_file_space_strategy(fcpl_id, H5F_FSPACE_STRATEGY_NONE, 0, 0);

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

    long number_of_cells = global_num_patches * num_cells_per_patch;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, vtkhdf_gid, "NumberOfCells", 1, dims, H5T_NATIVE_LONG, &number_of_cells);

    long number_of_points = global_num_patches * num_points_per_patch;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, vtkhdf_gid, "NumberOfPoints", 1, dims, H5T_NATIVE_LONG, &number_of_points);


    long number_of_connectivity_ids = global_num_patches * num_cells_per_patch * num_points_per_cell;
    dims[0] = 1;
    make_single_value_dataset_numerical(glob->mpirank, vtkhdf_gid, "NumberOfConnectivityIds", 1, dims, H5T_NATIVE_LONG, &number_of_connectivity_ids);


    int fits32 = number_of_points <= INT32_MAX
        && number_of_connectivity_ids <= INT32_MAX;

    hsize_t patch_dims[4] = {0,0,0,0};
    patch_dims[0] = num_cells_per_patch;
    make_dataset_numerical(glob, 
                           vtkhdf_gid, 
                           "Types", 
                           patch_dim, num_cells_per_patch_subdims,
                           1, patch_dims, 
                           num_patches_to_buffer, 
                           H5T_NATIVE_UINT8, 
                           types_cb,
                           &default_vtable);

    // write offsets
    patch_dims[0] = num_cells_per_patch;
    make_dataset_numerical(glob, 
                           vtkhdf_gid, 
                           "Offsets", 
                           patch_dim, num_cells_per_patch_subdims,
                           1, patch_dims, 
                           num_patches_to_buffer, 
                           fits32 ? H5T_NATIVE_INT32 : H5T_NATIVE_INT64,
                           fits32 ? get_offsets_int32_t : get_offsets_int64_t,
                           &offset_vtable);


    patch_dims[0] = num_points_per_patch;
    patch_dims[1] = 3;
    make_dataset_numerical(glob, 
                           vtkhdf_gid, 
                           "Points", 
                           patch_dim, num_points_per_patch_subdims,
                           2, patch_dims, 
                           num_patches_to_buffer, 
                           H5T_NATIVE_DOUBLE, 
                           coordinate_cb,
                           &default_vtable);

    patch_dims[0] = 1;
    if(patch_dim == 2)
    {
        patch_dims[1] = my+1;
        patch_dims[2] = mx+1;
    }
    else
    {
        patch_dims[1] = mz+1;
        patch_dims[2] = my+1;
        patch_dims[3] = mx+1;
    }
    make_dataset_numerical(glob, 
                           vtkhdf_gid, 
                           "PointIndexes", 
                           0, NULL,
                           patch_dim+1, patch_dims, 
                           num_patches_to_buffer, 
                           fits32 ? H5T_NATIVE_INT32 : H5T_NATIVE_INT64, 
                           fits32 ? get_point_indexes_int32_t : get_point_indexes_int64_t,
                           &default_vtable);

    patch_dims[0] = glob->domain->global_num_patches;
    hid_t src_space = H5Screate_simple(patch_dim+1, patch_dims, NULL);
    patch_dims[0] = glob->domain->global_num_patches * num_cells_per_patch * num_points_per_cell;
    hid_t virt_space = H5Screate_simple(1, patch_dims, NULL);

    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t slab_starts_3d[8][4] = {{0,0,0,0},{0,0,0,1},{0,0,1,1},{0,0,1,0},
                                    {0,1,0,0},{0,1,0,1},{0,1,1,1},{0,1,1,0}};
    hsize_t slab_starts_2d[4][3] = {{0,0,0},{0,0,1},{0,1,1},{0,1,0}};
    hsize_t block_dims_3d[4] = {1,mz,my,mx};
    hsize_t block_dims_2d[3] = {1,my,mx};
    for(hsize_t i=0; i < num_points_per_cell; i++)
    {
        hsize_t block_count[4] = {global_num_patches,1,1,1};
        status |= H5Sselect_hyperslab(src_space, 
                                      H5S_SELECT_SET, 
                                      patch_dim == 2 ? slab_starts_2d[i] : slab_starts_3d[i],
                                      NULL, 
                                      block_count, 
                                      patch_dim == 2 ? block_dims_2d : block_dims_3d);

        hsize_t stride = num_points_per_cell;
        hsize_t start = i;
        hsize_t count = number_of_cells;
        status |= H5Sselect_hyperslab(virt_space, H5S_SELECT_SET, &start, &stride, &count, NULL);
        status |= H5Pset_virtual(dcpl_id, virt_space, ".", "/VTKHDF/PointIndexes", src_space);
    }

    hid_t virt_id = H5Dcreate2(vtkhdf_gid, 
                               "Connectivity", 
                               fits32 ? H5T_NATIVE_INT32 : H5T_NATIVE_INT64,
                               virt_space, 
                               H5P_DEFAULT, 
                               dcpl_id, 
                               H5P_DEFAULT);
    status |= H5Dclose(virt_id);
    status |= H5Sclose(src_space);
    status |= H5Sclose(virt_space);

    /* avoid resource leaks by closing */
    status |= H5Gclose(vtkhdf_gid);


    hid_t celldata_gid = H5Gcreate2(file_id, celldata, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    patch_dims[0] = num_cells_per_patch;
    patch_dims[1] = clawpatch_opt->meqn;
    make_dataset_numerical(glob, 
                           celldata_gid, 
                           "meqn", 
                           patch_dim, num_cells_per_patch_subdims,
                           2, patch_dims, 
                           num_patches_to_buffer, 
                           H5T_NATIVE_FLOAT, 
                           value_cb,
                           &default_vtable);

    // write blockno
    patch_dims[0] = num_cells_per_patch;
    make_dataset_numerical(glob, 
                           celldata_gid, 
                           "blockno", 
                           patch_dim, num_cells_per_patch_subdims,
                           1, patch_dims, 
                           num_patches_to_buffer, 
                           H5T_NATIVE_INT, 
                           blockno_cb,
                           &default_vtable);

    //write patchno
    patch_dims[0] = num_cells_per_patch;
    make_dataset_numerical(glob, 
                           celldata_gid, 
                           "patchno", 
                           patch_dim, num_cells_per_patch_subdims,
                           1, patch_dims, 
                           num_patches_to_buffer, 
                           H5T_NATIVE_INT, 
                           patchno_cb,
                           &default_vtable);

    //write mpirank
    patch_dims[0] = num_cells_per_patch;
    make_dataset_numerical(glob, 
                           celldata_gid, 
                           "mpirank", 
                           patch_dim, num_cells_per_patch_subdims,
                           1, patch_dims, 
                           num_patches_to_buffer, 
                           H5T_NATIVE_INT, 
                           mpirank_cb,
                           &default_vtable);

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
                          coordinate_cb == NULL ? get_coordinates : coordinate_cb, 
                          coordinate_cb == NULL ? get_data : value_cb);
}

void fclaw_clawpatch_output_hdf5 (fclaw_global_t * glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d.vtkhdf", fclaw_opt->prefix, iframe);

    fclaw_clawpatch_output_hdf5_to_file(glob,basename, NULL, NULL);
}



