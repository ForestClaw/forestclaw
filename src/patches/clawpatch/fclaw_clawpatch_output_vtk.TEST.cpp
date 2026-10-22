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

#include <fclaw_clawpatch_output_vtk.h>
#include <test.hpp>

#include "fclaw_clawpatch_output_test_utils.hpp"

#include <cstring>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace
{

std::string
read_file(const std::string& filename)
{
    std::ifstream file(filename, std::ios::binary);
    if (!file)
    {
        return {};
    }
    return std::string(std::istreambuf_iterator<char>(file),
                       std::istreambuf_iterator<char>());
}

bool
vtu_dataarray_floats(const std::string& contents, const std::string& name,
                     std::vector<float>* values)
{
    const std::string name_attribute = "Name=\"" + name + "\"";
    const size_t tag_start = contents.find(name_attribute);
    if (tag_start == std::string::npos)
    {
        return false;
    }
    const size_t tag_end = contents.find('>', tag_start);
    const size_t offset_start = contents.find("offset=\"", tag_start);
    if (tag_end == std::string::npos || offset_start == std::string::npos ||
        offset_start > tag_end)
    {
        return false;
    }

    const size_t value_start = offset_start + std::strlen("offset=\"");
    const size_t value_end = contents.find('"', value_start);
    if (value_end == std::string::npos)
    {
        return false;
    }
    const size_t offset = std::stoull(contents.substr(value_start,
                                                     value_end - value_start));

    const std::string appended_start = "<AppendedData encoding=\"raw\">";
    const size_t appended_tag = contents.find(appended_start);
    if (appended_tag == std::string::npos)
    {
        return false;
    }
    const size_t underscore = contents.find('_', appended_tag + appended_start.size());
    if (underscore == std::string::npos)
    {
        return false;
    }

    const size_t data_start = underscore + 1 + offset;
    if (data_start + sizeof(uint64_t) > contents.size())
    {
        return false;
    }
    uint64_t byte_count;
    std::memcpy(&byte_count, contents.data() + data_start, sizeof(byte_count));
    if (byte_count % sizeof(float) != 0 ||
        data_start + sizeof(byte_count) + byte_count > contents.size())
    {
        return false;
    }
    values->resize(byte_count / sizeof(float));
    std::memcpy(values->data(), contents.data() + data_start + sizeof(byte_count),
                byte_count);
    return true;
}

} // namespace

TEST_CASE("VTU output contains the expected mesh and cell values")
{
    fclaw_clawpatch_output_test::TestData test_data;
    fclaw_clawpatch_output_vtk_to_file(test_data.glob, test_data.basename.c_str());

    const std::string contents = read_file(test_data.basename + ".vtu");
    REQUIRE_FALSE(contents.empty());
    CHECK(contents.find("<VTKFile type=\"UnstructuredGrid\"") != std::string::npos);
    CHECK(contents.find("<Piece NumberOfPoints=\"9\" NumberOfCells=\"4\">") !=
          std::string::npos);
    CHECK(contents.find("Name=\"meqn\"") != std::string::npos);

    std::vector<float> cell_values;
    REQUIRE(vtu_dataarray_floats(contents, "meqn", &cell_values));
    REQUIRE_EQ(cell_values.size(), 4);
    for (float value : cell_values)
    {
        CHECK_EQ(value, fclaw_clawpatch_output_test::output_value);
    }
}
