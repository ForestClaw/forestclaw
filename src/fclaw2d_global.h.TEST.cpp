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

#include <fclaw2d_global.h>
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

#ifdef FCLAW_ENABLE_DEBUG
    // TEST inserting twice
    CHECK_SC_ABORTED(fclaw_global_options_store(glob, key2, &option2));
#endif
    // Test with a non-existing key
    const char* key3 = "non-existing key";
#ifdef FCLAW_ENABLE_DEBUG
    CHECK_SC_ABORTED(fclaw_global_get_options(glob, key3));
#else
    void* retrieved_option3 = fclaw_global_get_options(glob, key3);
    CHECK_EQ(retrieved_option3, nullptr);
#endif

    fclaw_global_destroy(glob);
}

static bool destroyed;
static bool destroyed_2;
TEST_CASE("fclaw_global_attribute_store and fclaw2d_global_get_attribute test") {
    destroyed = false;
    destroyed_2 = false;
    fclaw_global_t* glob = fclaw_global_new();

    // Test with an integer
    int option1 = 10;
    const char* key1 = "option1";

    fclaw_global_attribute_store(glob, key1, &option1, [](void* data) {
        destroyed = true;
    });

    int* retrieved_option1 = static_cast<int*>(fclaw2d_global_get_attribute(glob, key1));
    CHECK_EQ(*retrieved_option1, option1);

    // Test with a string
    const char* option2 = "Test string";
    const char* key2 = "option2";
    fclaw_global_attribute_store(glob, key2, &option2, [](void* data) {
        destroyed_2 = true;
    });

    const char** retrieved_option2 = static_cast<const char**>(fclaw2d_global_get_attribute(glob, key2));
    CHECK_EQ(retrieved_option2, &option2);

#ifdef FCLAW_ENABLE_DEBUG
    // TEST inserting twice
    CHECK_SC_ABORTED(fclaw_global_attribute_store(glob, key2, &option2,nullptr));
#endif
    // Test with a non-existing key
    const char* key3 = "non-existing key";
#ifdef FCLAW_ENABLE_DEBUG
    CHECK_SC_ABORTED(fclaw2d_global_get_attribute(glob, key3));
#else
    void* retrieved_option3 = fclaw2d_global_get_attribute(glob, key3);
    CHECK_EQ(retrieved_option3, nullptr);
#endif

    fclaw_global_destroy(glob);
    CHECK_UNARY(destroyed);
    CHECK_UNARY(destroyed_2);
}

TEST_CASE("fclaw2d_global_set_global")
{
    fclaw_global_t* glob = (fclaw_global_t*)123;
    fclaw2d_global_set_global(glob);
    CHECK_EQ(fclaw2d_global_get_global(), glob);
    fclaw2d_global_unset_global();
}

TEST_CASE("fclaw2d_global_unset_global")
{
    fclaw_global_t* glob = (fclaw_global_t*)123;
    fclaw2d_global_set_global(glob);
    fclaw2d_global_unset_global();
#ifdef FCLAW_ENABLE_DEBUG
    CHECK_SC_ABORTED(fclaw2d_global_get_global());
#else
    CHECK_EQ(fclaw2d_global_get_global(), nullptr);
#endif
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fclaw2d_global_set_global twice fails")
{
    fclaw_global_t* glob = (fclaw_global_t*)123;
    fclaw2d_global_set_global(glob);
    CHECK_SC_ABORTED(fclaw2d_global_set_global(glob));
    fclaw2d_global_unset_global();
}

TEST_CASE("fclaw2d_global_unset_global assert fails when NULL")
{
    CHECK_SC_ABORTED(fclaw2d_global_unset_global());
}

TEST_CASE("fclaw2d_global_get_global assert fails when NULL")
{
    CHECK_SC_ABORTED(fclaw2d_global_get_global());
}

#endif
