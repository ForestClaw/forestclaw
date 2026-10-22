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

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include <fclaw_config.h>
#include <test.hpp>

#ifdef FCLAW_WITH_HDF5
#include <fclaw_clawpatch_output_hdf5.h>
#include <hdf5.h>

#include "fclaw_clawpatch_output_test_utils.hpp"

#include <cstdint>
#include <vector>

namespace
{

hsize_t
dataset_dimension(hid_t file, const char* name, int expected_rank, int dimension)
{
    hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
    REQUIRE(dataset >= 0);
    hid_t dataspace = H5Dget_space(dataset);
    REQUIRE(dataspace >= 0);

    const int rank = H5Sget_simple_extent_ndims(dataspace);
    REQUIRE_EQ(rank, expected_rank);
    std::vector<hsize_t> dimensions(static_cast<size_t>(rank));
    REQUIRE_EQ(H5Sget_simple_extent_dims(dataspace, dimensions.data(), nullptr),
               rank);

    H5Sclose(dataspace);
    H5Dclose(dataset);
    return dimensions[static_cast<size_t>(dimension)];
}

hsize_t
dataset_chunk_dimension(hid_t file, const char* name, int dimension)
{
    hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
    REQUIRE(dataset >= 0);
    hid_t properties = H5Dget_create_plist(dataset);
    REQUIRE(properties >= 0);

    hid_t dataspace = H5Dget_space(dataset);
    REQUIRE(dataspace >= 0);
    const int rank = H5Sget_simple_extent_ndims(dataspace);
    REQUIRE(rank > dimension);
    std::vector<hsize_t> dimensions(static_cast<size_t>(rank));
    REQUIRE_EQ(H5Pget_chunk(properties, rank, dimensions.data()), rank);

    H5Sclose(dataspace);
    H5Pclose(properties);
    H5Dclose(dataset);
    return dimensions[static_cast<size_t>(dimension)];
}

int64_t
read_hdf5_integer(hid_t file, const char* name)
{
    hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
    REQUIRE(dataset >= 0);
    int64_t value = -1;
    REQUIRE(H5Dread(dataset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &value) >= 0);
    H5Dclose(dataset);
    return value;
}

void
check_hdf5_values(hid_t file, const char* name, size_t count)
{
    hid_t dataset = H5Dopen2(file, name, H5P_DEFAULT);
    REQUIRE(dataset >= 0);
    std::vector<float> values(count);
    REQUIRE(H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    values.data()) >= 0);
    for (float value : values)
    {
        CHECK_EQ(value, fclaw_clawpatch_output_test::output_value);
    }
    H5Dclose(dataset);
}

} // namespace

TEST_CASE("VTKHDF output contains the expected mesh and cell values")
{
    fclaw_clawpatch_output_test::TestData test_data(sc_MPI_COMM_WORLD);
    test_data.opts->hdf5_compression_level = 5;
    const std::string filename = test_data.basename + ".vtkhdf";
    fclaw_clawpatch_output_hdf5_to_file(test_data.glob, filename.c_str(), nullptr,
                                        nullptr);

    unsigned long long local_patch_count = 0;
    for (int blockno = 0; blockno < test_data.domain->num_blocks; ++blockno)
    {
        local_patch_count += static_cast<unsigned long long>(
            test_data.domain->blocks[blockno].num_patches);
    }
    int mpisize = 1;
    REQUIRE(sc_MPI_Comm_size(test_data.comm, &mpisize) == 0);
    std::vector<unsigned long long> rank_patch_counts(
        static_cast<size_t>(mpisize));
    REQUIRE(sc_MPI_Allgather(&local_patch_count, 1, sc_MPI_UNSIGNED_LONG_LONG,
                             rank_patch_counts.data(), 1,
                             sc_MPI_UNSIGNED_LONG_LONG, test_data.comm) == 0);

    unsigned long long global_patch_count = 0;
    std::vector<unsigned long long> rank_starts(static_cast<size_t>(mpisize));
    for (int rank = 0; rank < mpisize; ++rank)
    {
        rank_starts[static_cast<size_t>(rank)] = global_patch_count * 9;
        global_patch_count += rank_patch_counts[static_cast<size_t>(rank)];
    }
    const hsize_t expected_points = static_cast<hsize_t>(global_patch_count * 9);
    const hsize_t expected_cells = static_cast<hsize_t>(global_patch_count * 4);

    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file >= 0);
    CHECK(read_hdf5_integer(file, "/VTKHDF/NumberOfPoints") == expected_points);
    CHECK(read_hdf5_integer(file, "/VTKHDF/NumberOfCells") == expected_cells);
    CHECK(dataset_dimension(file, "/VTKHDF/Points", 2, 0) == expected_points);
    CHECK(dataset_dimension(file, "/VTKHDF/Points", 2, 1) == 3);
    CHECK(dataset_dimension(file, "/VTKHDF/Connectivity", 1, 0) ==
          global_patch_count * 16);
    CHECK(dataset_dimension(file, "/VTKHDF/Offsets", 1, 0) ==
          global_patch_count * 4 + 1);
    CHECK(dataset_dimension(file, "/VTKHDF/CellData/meqn", 1, 0) ==
          expected_cells);
    const hsize_t points_chunk = dataset_chunk_dimension(file, "/VTKHDF/Points", 0);
    const hsize_t cells_chunk =
        dataset_chunk_dimension(file, "/VTKHDF/CellData/meqn", 0);
    const hsize_t offsets_chunk =
        dataset_chunk_dimension(file, "/VTKHDF/Offsets", 0);
    CHECK(points_chunk == 9);
    CHECK(cells_chunk == 4);
    CHECK(offsets_chunk == 4);
    for (unsigned long long start : rank_starts)
    {
        CHECK(start % points_chunk == 0);
        CHECK((start / 9 * 4) % cells_chunk == 0);
        CHECK((start / 9 * 4) % offsets_chunk == 0);
    }
    check_hdf5_values(file, "/VTKHDF/CellData/meqn",
                      static_cast<size_t>(expected_cells));
    CHECK(H5Lexists(file, "/VTKHDF/CellData/vtkGhostType", H5P_DEFAULT) == 0);
    H5Fclose(file);
}
#endif
