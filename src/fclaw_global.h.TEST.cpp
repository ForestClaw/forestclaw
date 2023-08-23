/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
#include <test.hpp>

TEST_CASE("fclaw_global_options_store and fclaw_global_get_options test") {
    fclaw_global_t* glob = fclaw_global_new();

    // Test with an integer
    int option1 = 10;
    const char* key1 = "option1";
    fclaw_global_options_store(glob, key1, &option1);

    int* retrieved_option1 = static_cast<int*>(fclaw_global_get_options(glob, key1));
    CHECK_EQ(*retrieved_option1, option1);

    // Test with a string
    const char* option2 = "Test string";
    const char* key2 = "option2";
    fclaw_global_options_store(glob, key2, &option2);

    const char** retrieved_option2 = static_cast<const char**>(fclaw_global_get_options(glob, key2));
    CHECK_EQ(retrieved_option2, &option2);

    // TEST inserting twice
    CHECK_SC_ABORTED(fclaw_global_options_store(glob, key2, &option2));
    const char* key3 = "non-existing key";
    CHECK_SC_ABORTED(fclaw_global_get_options(glob, key3));

    fclaw_global_destroy(glob);
}

static bool destroyed;
static bool destroyed_2;
TEST_CASE("fclaw_global_attribute_store and fclaw_global_get_attribute test") {
    destroyed = false;
    destroyed_2 = false;
    fclaw_global_t* glob = fclaw_global_new();

    // Test with an integer
    int attribute1 = 10;
    const char* key1 = "attribute1";

    fclaw_global_attribute_store(glob, key1, &attribute1, [](void* data) {
        destroyed = true;
    });

    int* retrieved_attribute1 = static_cast<int*>(fclaw_global_get_attribute(glob, key1));
    CHECK_EQ(*retrieved_attribute1, attribute1);

    // Test with a string
    const char* attribute2 = "Test string";
    const char* key2 = "attribute2";
    fclaw_global_attribute_store(glob, key2, &attribute2, [](void* data) {
        destroyed_2 = true;
    });

    const char** retrieved_attribute2 = static_cast<const char**>(fclaw_global_get_attribute(glob, key2));
    CHECK_EQ(retrieved_attribute2, &attribute2);

    fclaw_global_destroy(glob);
    CHECK_UNARY(destroyed);
    CHECK_UNARY(destroyed_2);
}

TEST_CASE("fclaw_global_set_static")
{
    fclaw_global_t* glob = (fclaw_global_t*)123;
    fclaw_global_set_static(glob);
    CHECK_EQ(fclaw_global_get_static_global(), glob);
    fclaw_global_clear_static();
}

TEST_CASE("fclaw_global_clear_static")
{
    fclaw_global_t* glob = (fclaw_global_t*)123;
    fclaw_global_set_static(glob);
    fclaw_global_clear_static();
    CHECK_EQ(fclaw_global_get_static_global(), nullptr);
}