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

#include <fclaw_pointer_map.h>
#include <fclaw_packing.h>
#include <fclaw_base.h>

#include <fclaw_global.h>
#include <test.hpp>

#include <initializer_list>
#include <vector>
TEST_CASE("fclaw_global_new default options")
{
    fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;


	CHECK_EQ(glob->curr_time, 0.0);
	CHECK_EQ(glob->curr_dt, 0.0);

	CHECK_EQ(glob->count_amr_advance, 0);
	CHECK_EQ(glob->count_ghost_exchange, 0);
	CHECK_EQ(glob->count_amr_regrid, 0);
	CHECK_EQ(glob->count_amr_new_domain, 0);
	CHECK_EQ(glob->count_single_step, 0);
	CHECK_EQ(glob->count_elliptic_grids, 0);
	CHECK_EQ(glob->count_multiproc_corner, 0);
	CHECK_EQ(glob->count_grids_per_proc, 0);
	CHECK_EQ(glob->count_grids_remote_boundary, 0);
	CHECK_EQ(glob->count_grids_local_boundary, 0);

	CHECK_EQ(glob->mpisize, 1);
	CHECK_EQ(glob->mpirank, 0);

	CHECK_EQ(fclaw_pointer_map_size(glob->vtables), 0);
	CHECK_EQ(fclaw_pointer_map_size(glob->options), 0);
	CHECK_EQ(fclaw_pointer_map_size(glob->attributes), 0);
	CHECK(glob->cont == nullptr);
	CHECK(glob->domain == nullptr);
	CHECK(glob->user == nullptr);


	fclaw_global_destroy(glob);
}

namespace
{

struct dummy_attribute
{
	size_t size;
	char value;

	dummy_attribute(size_t size_in, char value_in)
		:size(size_in),value(value_in)
	{
		//nothing to do
	}
};

size_t pack_dummy_options(fclaw_global_t *glob, void *user, char *buffer)
{
	dummy_attribute* options = (dummy_attribute*) user;

	char* buffer_start = buffer;
	buffer += fclaw_pack_size_t(options->size, buffer);
	for(size_t i = 0; i < options->size; i++){
		buffer[i] = options->value;
	}
	buffer += options->size;

	return buffer-buffer_start;
};
size_t unpack_dummy_options(fclaw_global_t *glob, char* buffer, void* user)
{
	dummy_attribute* attribute = (dummy_attribute*) user;

	char* buffer_start = buffer;

	size_t size;
	buffer += fclaw_unpack_size_t(buffer, &size);
	char value = buffer[0];
	for(size_t i = 1; i < size; i++){
		CHECK_EQ(buffer[i],value);
	}
	buffer += size;

	attribute->size = size;
	attribute->value = value;

	return buffer-buffer_start;
};
size_t packsize_dummy_attribute(fclaw_global_t *glob, void* user)
{
	dummy_attribute* options = (dummy_attribute*) user;
	return sizeof(size_t) + options->size;
};
void * new_dummy_attribute(fclaw_global_t *glob){
	return new dummy_attribute(0,0);
}
void destroy_dummy_attribute(void* user){
	delete (dummy_attribute*) user;
}

fclaw_packing_vtable_t dummy_opts_vt =
{
	pack_dummy_options,
	unpack_dummy_options,
	packsize_dummy_attribute,
	new_dummy_attribute,
	destroy_dummy_attribute
};

}
TEST_CASE("fclaw_global_pack with no packed attributes")
{
	for(bool extra_non_packed_attribute : {true, false})
	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	{
		fclaw_global_t* glob1;
    	glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;

		if(extra_non_packed_attribute)
		{
			dummy_attribute* attribute = new dummy_attribute(40, 'b');
			fclaw_global_attribute_store(glob1, "dummy_not_packed", attribute, NULL, destroy_dummy_attribute);
		}

		size_t packsize = fclaw_global_packsize(glob1);
		REQUIRE_GT(packsize, 0);

		char buffer[packsize];

		size_t bytes_written = fclaw_global_pack(glob1, buffer);

		REQUIRE_EQ(bytes_written, packsize);

		fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
		size_t bytes_read = fclaw_global_unpack(buffer, glob2);

		REQUIRE_EQ(bytes_read, packsize);

		CHECK_EQ(glob1->curr_time, glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,   glob2->curr_dt);

		REQUIRE_EQ(fclaw_pointer_map_size(glob2->attributes), 0);

		fclaw_global_destroy(glob1);
		fclaw_global_destroy(glob2);
	}
}
TEST_CASE("fclaw_global_pack with a single packed attribute")
{

	for(bool extra_non_packed_attribute : {true, false})
	for(bool already_stored_attribute : {true, false})
	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	{
		char dummy[] = "dummy1";
		std::vector<char*> args = {dummy};
		//fclaw_app_t* app = fclaw_app_new_on_comm(sc_MPI_COMM_WORLD,&argc,&argv,NULL);
		fclaw_global_t* glob1;
    	glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
		fclaw_global_vtable_store(glob1, "dummy1_vtable",  &dummy_opts_vt, NULL);

		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;

		dummy_attribute* attribute1 = new dummy_attribute(20, 'a');
		fclaw_global_attribute_store(glob1, "dummy1", attribute1, "dummy1_vtable", destroy_dummy_attribute);

		if(extra_non_packed_attribute)
		{
			dummy_attribute* attribute = new dummy_attribute(40, 'b');
			fclaw_global_attribute_store(glob1, "dummy_not_packed", attribute, NULL, destroy_dummy_attribute);
		}

		size_t packsize = fclaw_global_packsize(glob1);
		REQUIRE_GT(packsize, 0);

		char buffer[packsize];

		size_t bytes_written = fclaw_global_pack(glob1, buffer);

		REQUIRE_EQ(bytes_written, packsize);

		fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
		fclaw_global_vtable_store(glob2, "dummy1_vtable",  &dummy_opts_vt, NULL);

		dummy_attribute* attribute2 = NULL;
		if(already_stored_attribute)
		{
			attribute2 = new dummy_attribute(99, 'z');
			fclaw_global_attribute_store(glob2, "dummy1", attribute2, "dummy1_vtable", destroy_dummy_attribute);
		}

		size_t bytes_read = fclaw_global_unpack(buffer, glob2);

		attribute2 = (dummy_attribute*)(fclaw_global_get_attribute(glob2, "dummy1"));
		REQUIRE_NE(attribute2, nullptr);

		REQUIRE_EQ(bytes_read, packsize);

		CHECK_EQ(glob1->curr_time, glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,   glob2->curr_dt);

		REQUIRE_EQ(fclaw_pointer_map_size(glob2->attributes), 1);

		CHECK_EQ(attribute2->size, attribute1->size);
		CHECK_EQ(attribute2->value, attribute1->value);

		fclaw_global_destroy(glob1);
		fclaw_global_destroy(glob2);
	}
}
TEST_CASE("fclaw_global_pack with two packed attributes")
{

	for(bool extra_non_packed_attribute : {true, false})
	for(bool already_stored_attribute_1 : {true, false})
	for(bool already_stored_attribute_2 : {true, false})
	for(double curr_time : {1.0, 1.2})
	for(double curr_dt : {1.0, 1.2})
	{
		fclaw_global_t* glob1;
    	glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
		fclaw_global_vtable_store(glob1, "dummy1",  &dummy_opts_vt, NULL);
		fclaw_global_vtable_store(glob1, "dummy2",  &dummy_opts_vt, NULL);
		glob1->curr_time                    = curr_time;
		glob1->curr_dt                      = curr_dt;

		dummy_attribute* attribute1_1 = new dummy_attribute(20, 'a');
		dummy_attribute* attribute1_2 = new dummy_attribute(40, 'b');

		fclaw_global_attribute_store(glob1, "dummy1", attribute1_1, "dummy1", destroy_dummy_attribute);
		fclaw_global_attribute_store(glob1, "dummy2", attribute1_2, "dummy2", destroy_dummy_attribute);

		if(extra_non_packed_attribute)
		{
			dummy_attribute* attribute = new dummy_attribute(40, 'b');
			fclaw_global_attribute_store(glob1, "dummy_not_packed", attribute, NULL, destroy_dummy_attribute);
		}

		size_t packsize = fclaw_global_packsize(glob1);
		REQUIRE_GT(packsize, 0);

		char buffer[packsize];

		size_t bytes_written = fclaw_global_pack(glob1, buffer);

		REQUIRE_EQ(bytes_written, packsize);

		fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
		fclaw_global_vtable_store(glob2, "dummy1",  &dummy_opts_vt, NULL);
		fclaw_global_vtable_store(glob2, "dummy2",  &dummy_opts_vt, NULL);

		dummy_attribute* attribute2_1 = NULL;
		dummy_attribute* attribute2_2 = NULL;
		if(already_stored_attribute_1)
		{
			attribute2_1 = new dummy_attribute(99, 'z');
			fclaw_global_attribute_store(glob2, "dummy1", attribute2_1, "dummy1", destroy_dummy_attribute);
		}
		if(already_stored_attribute_2)
		{
			attribute2_2 = new dummy_attribute(99, 'z');
			fclaw_global_attribute_store(glob2, "dummy2", attribute2_2, "dummy2", destroy_dummy_attribute);
		}

		size_t bytes_read = fclaw_global_unpack(buffer, glob2);

		attribute2_1 = (dummy_attribute*)(fclaw_global_get_attribute(glob2, "dummy1"));
		REQUIRE_NE(attribute2_1, nullptr);
		attribute2_2 = (dummy_attribute*)(fclaw_global_get_attribute(glob2, "dummy2"));
		REQUIRE_NE(attribute2_2, nullptr);

		REQUIRE_EQ(bytes_read, packsize);

		CHECK_EQ(glob1->curr_time, glob2->curr_time);
		CHECK_EQ(glob1->curr_dt,   glob2->curr_dt);

		REQUIRE_EQ(fclaw_pointer_map_size(glob2->attributes), 2);

		CHECK_EQ(attribute2_1->size, attribute1_1->size);
		CHECK_EQ(attribute2_1->value, attribute1_1->value);

		CHECK_EQ(attribute2_2->size, attribute1_2->size);
		CHECK_EQ(attribute2_2->value, attribute1_2->value);

		fclaw_global_destroy(glob1);
		fclaw_global_destroy(glob2);
	}
}
TEST_CASE("fclaw_global_pack aborts with unregistered vtable")
{
	fclaw_global_t* glob1;
   	glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	glob1->curr_time                    = 100;
	glob1->curr_dt                      = 200;

	dummy_attribute* attribute = new dummy_attribute(20, 'a');
	fclaw_global_attribute_store(glob1, "dummy1", attribute, "pack_dummy1", destroy_dummy_attribute);

	char buffer[BUFSIZ];
	CHECK_SC_ABORTED(fclaw_global_pack(glob1, buffer));
}
TEST_CASE("fclaw_global_packsize aborts with unregistered vtable")
{
	fclaw_global_t* glob1;
   	glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	glob1->curr_time                    = 100;
	glob1->curr_dt                      = 200;

	dummy_attribute* attribute = new dummy_attribute(20, 'a');
	fclaw_global_attribute_store(glob1, "dummy1", attribute, "pack_dummy1", destroy_dummy_attribute);

	CHECK_SC_ABORTED(fclaw_global_packsize(glob1));
}
TEST_CASE("fclaw_global_unpack aborts with unregistered vtable")
{
	fclaw_global_t* glob1;
   	glob1 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	glob1->curr_time                    = 1;
	glob1->curr_dt                      = 1;

	dummy_attribute* attribute = new dummy_attribute(20, 'a');
	fclaw_global_vtable_store(glob1, "pack_dummy1",  &dummy_opts_vt, NULL);
	fclaw_global_attribute_store(glob1, "dummy1", attribute, "pack_dummy1", destroy_dummy_attribute);

	size_t packsize = fclaw_global_packsize(glob1);
	REQUIRE_GT(packsize, 0);

	char buffer[packsize];

	size_t bytes_written = fclaw_global_pack(glob1, buffer);

	REQUIRE_EQ(bytes_written, packsize);

	fclaw_global_t* glob2=fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;
	//no vtables in glob2
	CHECK_SC_ABORTED(fclaw_global_unpack(buffer, glob2));
}

TEST_CASE("fclaw_global_options_store and fclaw_global_get_options test") {
    fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;

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
	const char* option3 = "Test string 2";
	fclaw_global_options_store(glob, key2, &option3);
	CHECK_EQ(fclaw_global_get_options(glob, key2), &option3);

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

TEST_CASE("fclaw_global_set_global")
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
#ifdef FCLAW_ENABLE_DEBUG
    CHECK_SC_ABORTED(fclaw_global_get_static_global());
#else
    CHECK_EQ(fclaw_global_get_static_global(), nullptr);
#endif
}

#ifdef FCLAW_ENABLE_DEBUG

TEST_CASE("fclaw_global_set_global twice fails")
{
    fclaw_global_t* glob = (fclaw_global_t*)123;
    fclaw_global_set_static(glob);
    CHECK_SC_ABORTED(fclaw_global_set_static(glob));
    fclaw_global_clear_static();
}

TEST_CASE("fclaw_global_unset_global assert fails when NULL")
{
    CHECK_SC_ABORTED(fclaw_global_clear_static());
}

TEST_CASE("fclaw_global_get_global assert fails when NULL")
{
    CHECK_SC_ABORTED(fclaw_global_get_static_global());
}

#endif

TEST_CASE("fclaw_global_vtable_store and fclaw_global_get_vtable test") {
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;

	// Test with an integer
	int vtable1 = 10;
	const char* key1 = "vtable1";
	fclaw_global_vtable_store(glob, key1, &vtable1,  nullptr);

	int* retrieved_vtable1 = static_cast<int*>(fclaw_global_get_vtable(glob, key1));
	CHECK_EQ(*retrieved_vtable1, vtable1);

	// Test with a string
	const char* vtable2 = "Test string";
	const char* key2 = "vtable2";
	fclaw_global_vtable_store(glob, key2, &vtable2,  nullptr);

	const char** retrieved_vtable2 = static_cast<const char**>(fclaw_global_get_vtable(glob, key2));
	CHECK_EQ(*retrieved_vtable2, vtable2);

	// Test with a non-existing key
	const char* key3 = "non-existing key";
	void* retrieved_vtable3 = fclaw_global_get_vtable(glob, key3);
	CHECK_EQ(retrieved_vtable3, nullptr);

	// Test with a destroy callback
	int* vtable4 = new int(20);
	const char* key4 = "vtable4";
	fclaw_global_vtable_store(glob, key4, vtable4, [](void* ptr) { delete static_cast<int*>(ptr); });

	int* retrieved_vtable4 = static_cast<int*>(fclaw_global_get_vtable(glob, key4));
	CHECK_EQ(*retrieved_vtable4, *vtable4);

	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_global_attribute_store and fclaw_global_get_attribute test") {
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);;

	// Test with an integer
	int attribute1 = 10;
	const char* key1 = "attribute1";
	fclaw_global_attribute_store(glob, key1, &attribute1, nullptr, nullptr);

	int* retrieved_attribute1 = static_cast<int*>(fclaw_global_get_attribute(glob, key1));
	CHECK_EQ(*retrieved_attribute1, attribute1);

	// Test with a string
	const char* attribute2 = "Test string";
	const char* key2 = "attribute2";
	fclaw_global_attribute_store(glob, key2, &attribute2, nullptr, nullptr);

	const char** retrieved_attribute2 = static_cast<const char**>(fclaw_global_get_attribute(glob, key2));
	CHECK_EQ(*retrieved_attribute2, attribute2);

	// Test with a non-existing key
	const char* key3 = "non-existing key";
	void* retrieved_attribute3 = fclaw_global_get_attribute(glob, key3);
	CHECK_EQ(retrieved_attribute3, nullptr);

	// Test with a destroy callback
	int* attribute4 = new int(20);
	const char* key4 = "attribute4";
	fclaw_global_attribute_store(glob, key4, attribute4, nullptr, [](void* ptr) { delete static_cast<int*>(ptr); });

	int* retrieved_attribute4 = static_cast<int*>(fclaw_global_get_attribute(glob, key4));
	CHECK_EQ(*retrieved_attribute4, *attribute4);

	fclaw_global_destroy(glob);
}

