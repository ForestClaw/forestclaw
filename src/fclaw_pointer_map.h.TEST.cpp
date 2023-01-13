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

#include <fclaw_pointer_map.h>
#include <test.hpp>
#include <map>
#include <string>

TEST_CASE("fclaw_pointer_map empty map get is nullptr")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	CHECK_EQ(fclaw_pointer_map_get(map, "key"), nullptr);
	fclaw_pointer_map_destroy(map);
}

TEST_CASE("fclaw_pointer_map set and get null destory")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	void* value = (void*) 1;
	fclaw_pointer_map_insert(map, "key", value, nullptr);
	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value);
	fclaw_pointer_map_destroy(map);
}

bool value_destroyed = false;
TEST_CASE("fclaw_pointer_map set and get")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	void* value = (void*) 1;
	value_destroyed = false;
	fclaw_pointer_map_insert(map, "key", value, [](void* v) {
		value_destroyed = true;
	});
	CHECK_UNARY_FALSE(value_destroyed);
	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value);
	fclaw_pointer_map_destroy(map);
	CHECK_UNARY(value_destroyed);
}

TEST_CASE("fclaw_pointer_map set and get modified key")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	void* value = (void*) 1;
	value_destroyed = false;
	std::string key = "key";
	fclaw_pointer_map_insert(map, key.c_str(), value, [](void* v) {
		value_destroyed = true;
	});
	key = "blah";
	CHECK_UNARY_FALSE(value_destroyed);
	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value);
	fclaw_pointer_map_destroy(map);
	CHECK_UNARY(value_destroyed);
}

bool value2_destroyed = false;
TEST_CASE("fclaw_pointer_map set and get two values")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	
	void* value = (void*) 1;
	void* value2 = (void*) 2;

	value_destroyed = false;
	value2_destroyed = false;

	fclaw_pointer_map_insert(map, "key", value, [](void* v) {
		value_destroyed = true;
	});
	CHECK_UNARY_FALSE(value_destroyed);

	fclaw_pointer_map_insert(map, "key2", value2, [](void* v) {
		value2_destroyed = true;
	});
	CHECK_UNARY_FALSE(value2_destroyed);

	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value);
	CHECK_EQ(fclaw_pointer_map_get(map, "key2"), value2);

	fclaw_pointer_map_destroy(map);
	CHECK_UNARY(value_destroyed);
	CHECK_UNARY(value2_destroyed);
}

TEST_CASE("fclaw_pointer_map set twice and get")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	
	void* value = (void*) 1;
	void* value2 = (void*) 2;

	value_destroyed = false;
	value2_destroyed = false;

	fclaw_pointer_map_insert(map, "key", value, [](void* v) {
		value_destroyed = true;
	});
	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value);

	CHECK_UNARY_FALSE(value_destroyed);

	fclaw_pointer_map_insert(map, "key", value2, [](void* v) {
		value2_destroyed = true;
	});

	CHECK_UNARY_FALSE(value2_destroyed);
	//old value should be destroyed
	CHECK_UNARY(value_destroyed);

	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value2);

	fclaw_pointer_map_destroy(map);
	CHECK_UNARY(value2_destroyed);
}

TEST_CASE("fclaw_pointer_map set twice and get with null destory")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	
	void* value = (void*) 1;
	void* value2 = (void*) 2;

	fclaw_pointer_map_insert(map, "key", value, nullptr);

	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value);

	fclaw_pointer_map_insert(map, "key", value2, nullptr);

	CHECK_EQ(fclaw_pointer_map_get(map, "key"), value2);

	fclaw_pointer_map_destroy(map);
}

TEST_CASE("fclaw_pointer_map iterate on empty map")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	fclaw_pointer_map_iterate(map,
	                         [](const char* key, void* value, void* user){
								CHECK_UNARY(false);
	                          },
	                          NULL);
	fclaw_pointer_map_destroy(map);
}

TEST_CASE("fclaw_pointer_map iterate on map with two values")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	void* value = (void*) 1;
	void* value2 = (void*) 2;

	fclaw_pointer_map_insert(map, "a", value, nullptr);

	fclaw_pointer_map_insert(map, "b", value2, nullptr);


	std::map<std::string,void*> iterated_values;

	fclaw_pointer_map_iterate(map,
	                         [](const char* key, void* value, void* user){
								std::map<std::string,void*>& iterated_values = *(std::map<std::string,void*>*)user;
								CHECK_EQ(iterated_values.count(key),0);
								iterated_values[key] = value;
	                          },
	                          &iterated_values);
	CHECK_EQ(iterated_values["a"],value);
	CHECK_EQ(iterated_values["b"],value2);
	fclaw_pointer_map_destroy(map);
}

TEST_CASE("fclaw_pointer_size on empty map")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	CHECK_EQ(fclaw_pointer_map_size(map), 0);
	fclaw_pointer_map_destroy(map);
}

TEST_CASE("fclaw_pointer_size with two values")
{
	fclaw_pointer_map_t* map = fclaw_pointer_map_new();
	void* value = (void*) 1;
	void* value2 = (void*) 2;

	fclaw_pointer_map_insert(map, "a", value, nullptr);

	fclaw_pointer_map_insert(map, "b", value2, nullptr);


	CHECK_EQ(fclaw_pointer_map_size(map), 2);

	fclaw_pointer_map_destroy(map);
}