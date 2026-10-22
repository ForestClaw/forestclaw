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
    fclaw_clawpatch_output_test::TestData test_data;
    const std::string filename = test_data.basename + ".vtkhdf";
    fclaw_clawpatch_output_hdf5_to_file(test_data.glob, filename.c_str(), nullptr,
                                        nullptr);

    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    REQUIRE(file >= 0);
    CHECK(read_hdf5_integer(file, "/VTKHDF/NumberOfPoints") == 9);
    CHECK(read_hdf5_integer(file, "/VTKHDF/NumberOfCells") == 4);
    CHECK(dataset_dimension(file, "/VTKHDF/Points", 2, 0) == 9);
    CHECK(dataset_dimension(file, "/VTKHDF/Points", 2, 1) == 3);
    CHECK(dataset_dimension(file, "/VTKHDF/Connectivity", 1, 0) == 16);
    CHECK(dataset_dimension(file, "/VTKHDF/Offsets", 1, 0) == 5);
    CHECK(dataset_dimension(file, "/VTKHDF/CellData/meqn", 1, 0) == 4);
    check_hdf5_values(file, "/VTKHDF/CellData/meqn", 4);
    CHECK(H5Lexists(file, "/VTKHDF/CellData/vtkGhostType", H5P_DEFAULT) == 0);
    H5Fclose(file);
}
#endif
