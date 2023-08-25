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
#include <fclaw_packing.h>
#include <vector>
#include <fclaw2d_convenience.h>

namespace
{

struct dummy_options
{
	size_t size;
	char value;

	dummy_options(size_t size_in, char value_in)
		:size(size_in),value(value_in)
	{
		//nothing to do
	}
};

size_t pack_dummy_options(void* user, char* buffer)
{
	dummy_options* options = (dummy_options*) user;

	char* buffer_start = buffer;
	buffer += fclaw_pack_size_t(options->size, buffer);
	for(size_t i = 0; i < options->size; i++){
		buffer[i] = options->value;
	}
	buffer += options->size;

	return buffer-buffer_start;
};
size_t unpack_dummy_options(char* buffer, void** user)
{

	char* buffer_start = buffer;

	size_t size;
	buffer += fclaw_unpack_size_t(buffer, &size);
	char value = buffer[0];
	for(size_t i = 1; i < size; i++){
		CHECK_EQ(buffer[i],value);
	}
	buffer += size;

	*user = new dummy_options(size,value);

	return buffer-buffer_start;
};
size_t packsize_dummy_options(void* user)
{
	dummy_options* options = (dummy_options*) user;
	return sizeof(size_t) + options->size;
};

void destroy_dummy_options(void* user){
	delete (dummy_options*) user;
}

fclaw_packing_vtable_t dummy_opts_vt =
{
	pack_dummy_options,
	unpack_dummy_options,
	packsize_dummy_options,
	destroy_dummy_options
};

}
TEST_CASE("fclaw_global_pack with no options structure")
{

	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	{
		fclaw_global_t* glob1;
    	glob1 = fclaw_global_new();
		glob1->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;

		size_t packsize = fclaw_global_packsize(glob1);
		REQUIRE_GT(packsize, 0);

		char buffer[packsize];

		size_t bytes_written = fclaw_global_pack(glob1, buffer);

		REQUIRE_EQ(bytes_written, packsize);

		fclaw_global_t* glob2;
		size_t bytes_read = fclaw_global_unpack(buffer, &glob2);

		REQUIRE_EQ(bytes_read, packsize);

		CHECK_EQ(glob1->curr_time, glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,   glob2->curr_dt);

		REQUIRE_EQ(fclaw_pointer_map_size(glob2->options), 0);

		fclaw_global_destroy(glob1);
		fclaw_global_destroy(glob2);
	}
}
TEST_CASE("fclaw_global_pack with a single options structure")
{

	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	{
		char dummy[] = "dummy1";
		std::vector<char*> args = {dummy};
		char ** argv = args.data();
		int argc = 1;
		//fclaw_app_t* app = fclaw_app_new_on_comm(sc_MPI_COMM_WORLD,&argc,&argv,NULL);
		fclaw_global_t* glob1;
    	glob1 = fclaw_global_new();
		glob1->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;

		dummy_options* options = new dummy_options(20, 'a');
		fclaw_app_register_options_packing_vtable("dummy1",  &dummy_opts_vt);
		fclaw_pointer_map_insert(glob1->options, "dummy1", options, destroy_dummy_options);

		size_t packsize = fclaw_global_packsize(glob1);
		REQUIRE_GT(packsize, 0);

		char buffer[packsize];

		size_t bytes_written = fclaw_global_pack(glob1, buffer);

		REQUIRE_EQ(bytes_written, packsize);

		fclaw_global_t* glob2;
		size_t bytes_read = fclaw_global_unpack(buffer, &glob2);

		REQUIRE_EQ(bytes_read, packsize);

		CHECK_EQ(glob1->curr_time, glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,   glob2->curr_dt);

		REQUIRE_EQ(fclaw_pointer_map_size(glob2->options), 1);

		dummy_options* options_unpacked = (dummy_options*) fclaw_pointer_map_get(glob2->options, "dummy1");
		REQUIRE(options_unpacked != nullptr);

		CHECK_EQ(options_unpacked->size, options->size);
		CHECK_EQ(options_unpacked->value, options->value);



		fclaw_global_destroy(glob1);
		fclaw_global_destroy(glob2);
	}
}
TEST_CASE("fclaw_global_pack with two options structure")
{

	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	{
		fclaw_global_t* glob1;
    	glob1 = fclaw_global_new();
		glob1->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;

		dummy_options* options = new dummy_options(20, 'a');
		dummy_options* options2 = new dummy_options(40, 'b');
		fclaw_app_register_options_packing_vtable("dummy1",  &dummy_opts_vt);
		fclaw_pointer_map_insert(glob1->options, "dummy1", options, destroy_dummy_options);
		fclaw_app_register_options_packing_vtable("dummy2",  &dummy_opts_vt);
		fclaw_pointer_map_insert(glob1->options, "dummy2", options2, destroy_dummy_options);

		size_t packsize = fclaw_global_packsize(glob1);
		REQUIRE_GT(packsize, 0);

		char buffer[packsize];

		size_t bytes_written = fclaw_global_pack(glob1, buffer);

		REQUIRE_EQ(bytes_written, packsize);

		fclaw_global_t* glob2;
		size_t bytes_read = fclaw_global_unpack(buffer, &glob2);

		REQUIRE_EQ(bytes_read, packsize);

		CHECK_EQ(glob1->curr_time, glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,   glob2->curr_dt);

		REQUIRE_EQ(fclaw_pointer_map_size(glob2->options), 2);

		dummy_options* options_unpacked = (dummy_options*) fclaw_pointer_map_get(glob2->options, "dummy1");
		dummy_options* options_unpacked2 = (dummy_options*) fclaw_pointer_map_get(glob2->options, "dummy2");
		REQUIRE(options_unpacked != nullptr);
		REQUIRE(options_unpacked2 != nullptr);

		CHECK_EQ(options_unpacked->size, options->size);
		CHECK_EQ(options_unpacked->value, options->value);

		CHECK_EQ(options_unpacked2->size, options2->size);
		CHECK_EQ(options_unpacked2->value, options2->value);

		fclaw_global_destroy(glob1);
		fclaw_global_destroy(glob2);
	}
}
TEST_CASE("fclaw_global_pack aborts with unregistered vtable")
{
	fclaw_global_t* glob1;
   	glob1 = fclaw_global_new();
	glob1->curr_time                    = 100;
	glob1->curr_dt                      = 200;
	glob1->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);

	dummy_options* options = new dummy_options(20, 'a');
	fclaw_pointer_map_insert(glob1->options, "dummy1", options, destroy_dummy_options);

	char buffer[100];
	CHECK_SC_ABORTED(fclaw_global_pack(glob1, buffer));
}
TEST_CASE("fclaw_global_packsize aborts with unregistered vtable")
{
	fclaw_global_t* glob1;
   	glob1 = fclaw_global_new();
	glob1->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);
	glob1->curr_time                    = 100;
	glob1->curr_dt                      = 200;

	dummy_options* options = new dummy_options(20, 'a');
	fclaw_pointer_map_insert(glob1->options, "dummy1", options, destroy_dummy_options);

	CHECK_SC_ABORTED(fclaw_global_packsize(glob1));
}
TEST_CASE("fclaw_global_unppack aborts with unregistered vtable")
{
	fclaw_global_t* glob1;
   	glob1 = fclaw_global_new();
	glob1->curr_time                    = 1;
	glob1->curr_dt                      = 1;

	glob1->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);

	dummy_options* options = new dummy_options(20, 'a');
	fclaw_app_register_options_packing_vtable("dummy1",  &dummy_opts_vt);
	fclaw_pointer_map_insert(glob1->options, "dummy1", options, destroy_dummy_options);

	size_t packsize = fclaw_global_packsize(glob1);
	REQUIRE_GT(packsize, 0);

	char buffer[packsize];

	size_t bytes_written = fclaw_global_pack(glob1, buffer);

	REQUIRE_EQ(bytes_written, packsize);

	fclaw_global_t* glob2=nullptr;
	fclaw_app_register_options_packing_vtable("dummy1",  nullptr);
	CHECK_SC_ABORTED(fclaw_global_unpack(buffer, &glob2));
}
TEST_CASE("fclaw_global_options_store and fclaw_global_get_options test") {
    fclaw_global_t* glob = fclaw_global_new();
	glob->domain = fclaw2d_domain_new_unitsquare(sc_MPI_COMM_WORLD, 0);

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