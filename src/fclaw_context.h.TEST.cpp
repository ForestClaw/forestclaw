/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "fclaw_global.h"
#include <fclaw_context.h>
#include <test.hpp>

TEST_CASE("fclaw_context_get new context")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context = fclaw_context_get(glob, "test");
	CHECK_NE(context, nullptr);
	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_context_get two new contexts")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context1 = fclaw_context_get(glob, "test1");
	fclaw_context_t *context2 = fclaw_context_get(glob, "test2");
	CHECK_NE(context1, nullptr);
	CHECK_NE(context2, nullptr);
	CHECK_NE(context1, context2);
	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_context_get existing context without saving")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context1 = fclaw_context_get(glob, "test");
	CHECK_NE(context1, nullptr);
	//should fail since it wasnn't saved
	CHECK_SC_ABORTED(fclaw_context_get(glob, "test"));
}

TEST_CASE("fclaw_context_get existing context")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context1 = fclaw_context_get(glob, "test");
	CHECK_NE(context1, nullptr);
	fclaw_context_save(context1);
	fclaw_context_t *context2 = fclaw_context_get(glob, "test");
	CHECK_EQ(context1, context2);
	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_context_get two existing contexts")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context1 = fclaw_context_get(glob, "test1");
	fclaw_context_t *context2 = fclaw_context_get(glob, "test2");
	CHECK_NE(context1, nullptr);
	CHECK_NE(context2, nullptr);
	CHECK_NE(context1, context2);
	fclaw_context_save(context1);
	fclaw_context_save(context2);
	fclaw_context_t *context3 = fclaw_context_get(glob, "test1");
	fclaw_context_t *context4 = fclaw_context_get(glob, "test2");
	CHECK_EQ(context1, context3);
	CHECK_EQ(context2, context4);
	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_context_get_int new context")
{
	for(int default_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);
		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int called twice on new context")
{
	for(int default_value : {-100, 0, 42})
	for(int changed_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);

		value = changed_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, changed_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int new context two values")
{
	for(int default_value1 : {-100, 0, 42})
	for(int default_value2 : {-1, 1, 1024})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		int value1 = default_value1;
		fclaw_context_get_int(context, "test1", &value1);
		CHECK_EQ(value1, default_value1);
		int value2 = default_value2;
		fclaw_context_get_int(context, "test2", &value2);
		CHECK_EQ(value2, default_value2);
		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int existing context")
{
	for(int default_value : {-100, 0, 42})
	for(int changed_value : {-1, 4, 82})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value1 =  default_value;

		fclaw_context_get_int(context, "test", &value1);
		CHECK_EQ(value1, default_value);

		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");

		int value2 = 0;
		fclaw_context_get_int(context, "test", &value2);
		CHECK_EQ(value2, default_value);

		value2 = changed_value;
		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");

		int value3 = 0;
		fclaw_context_get_int(context, "test", &value3);
		CHECK_EQ(value3, changed_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int existing context two values")
{
	for(int default_value1 : {-100, 0, 42})
	for(int changed_value1 : {-1, 4, 82})
	for(int default_value2: {-100, 0, 42})
	for(int changed_value2: {-1, 4, 82})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value = default_value1;
		int value2 = default_value2;

		fclaw_context_get_int(context, "test1", &value);
		fclaw_context_get_int(context, "test2", &value2);

		CHECK_EQ(value, default_value1);
		CHECK_EQ(value2, default_value2);

		fclaw_context_save(context);

		fclaw_context_get_int(context, "test1", &value);
		fclaw_context_get_int(context, "test2", &value2);

		CHECK_EQ(value, default_value1);
		CHECK_EQ(value2, default_value2);

		value = changed_value1;
		value2 = changed_value2;

		fclaw_context_save(context);

		fclaw_context_get_int(context, "test1", &value);
		fclaw_context_get_int(context, "test2", &value2);

		CHECK_EQ(value, changed_value1);
		CHECK_EQ(value2, changed_value2);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int called for non-existing value")
{
	for(int default_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		fclaw_context_save(context);

		// second get call, since we didn't get_int in first, should fail
		context = fclaw_context_get(glob, "test");
		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int called for non-exising value other value")
{
	for(int default_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		int value = 0;
		fclaw_context_get_int(context, "test", &value);
		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");
		int value2 = default_value;
		fclaw_context_get_int(context, "test-does-not-exist", &value2);
		CHECK_EQ(value2, default_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_int save without getting all variables")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context = fclaw_context_get(glob, "test");
	int value = 0;
	fclaw_context_get_int(context, "test", &value);
	fclaw_context_save(context);

	// second get call
	context = fclaw_context_get(glob, "test");
	// change value, shouldn't save, since get wasn't called
	value = 1;
	fclaw_context_save(context);

	// check that value wasn't changed
	context = fclaw_context_get(glob, "test");
	fclaw_context_get_int(context, "test", &value);
	CHECK_EQ(value, 0);

	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_context_get_double new context")
{
	for(double default_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		double value = default_value;
		fclaw_context_get_double(context, "test", &value);
		CHECK_EQ(value, default_value);
		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double called twice on new context")
{
	for(double default_value : {-100, 0, 42})
	for(double changed_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		double value = default_value;
		fclaw_context_get_double(context, "test", &value);
		CHECK_EQ(value, default_value);

		value = changed_value;
		fclaw_context_get_double(context, "test", &value);
		CHECK_EQ(value, changed_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double new context two values")
{
	for(double default_value1 : {-100, 0, 42})
	for(double default_value2 : {-1, 1, 1024})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		double value1 = default_value1;
		fclaw_context_get_double(context, "test1", &value1);
		CHECK_EQ(value1, default_value1);
		double value2 = default_value2;
		fclaw_context_get_double(context, "test2", &value2);
		CHECK_EQ(value2, default_value2);
		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double existing context")
{
	for(double default_value : {-100, 0, 42})
	for(double changed_value : {-1, 4, 82})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		double value = default_value;

		fclaw_context_get_double(context, "test", &value);
		CHECK_EQ(value, default_value);

		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");

		double value2 = 0;
		fclaw_context_get_double(context, "test", &value2);
		CHECK_EQ(value2, default_value);

		value2 = changed_value;

		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");

		double value3 = 0;
		fclaw_context_get_double(context, "test", &value3);
		CHECK_EQ(value3, changed_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double existing context two values")
{
	for(double default_value1 : {-100, 0, 42})
	for(double changed_value1 : {-1, 4, 82})
	for(double default_value2: {-100, 0, 42})
	for(double changed_value2: {-1, 4, 82})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		double value = default_value1;
		double value2 = default_value2;

		fclaw_context_get_double(context, "test1", &value);
		fclaw_context_get_double(context, "test2", &value2);

		CHECK_EQ(value, default_value1);
		CHECK_EQ(value2, default_value2);

		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");

		value = 0; value2 = 0;

		fclaw_context_get_double(context, "test1", &value);
		fclaw_context_get_double(context, "test2", &value2);

		CHECK_EQ(value, default_value1);
		CHECK_EQ(value2, default_value2);

		value = changed_value1;
		value2 = changed_value2;

		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");

		value = 0; value2 = 0;

		fclaw_context_get_double(context, "test1", &value);
		fclaw_context_get_double(context, "test2", &value2);

		CHECK_EQ(value, changed_value1);
		CHECK_EQ(value2, changed_value2);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double called for non-existing value")
{
	for(double default_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		fclaw_context_save(context);

		// second get call, since we didn't get_double in first, should fail
		context = fclaw_context_get(glob, "test");
		double value = default_value;
		fclaw_context_get_double(context, "test", &value);

		CHECK_EQ(value, default_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double called for non-exising value other value")
{
	for(double default_value : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");
		double value = 0;
		fclaw_context_get_double(context, "test", &value);
		fclaw_context_save(context);

		context = fclaw_context_get(glob, "test");
		double value2 = default_value;
		fclaw_context_get_double(context, "test-does-not-exist", &value2);
		CHECK_EQ(value2, default_value);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double save without getting all variables")
{
	fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
	fclaw_context_t *context = fclaw_context_get(glob, "test");
	double value = 0;
	fclaw_context_get_double(context, "test", &value);
	fclaw_context_save(context);

	// second get call
	context = fclaw_context_get(glob, "test");
	// change value, shouldn't save, since get wasn't called
	value = 1;
	fclaw_context_save(context);
	
	// check that value wasn't changed
	context = fclaw_context_get(glob, "test");
	fclaw_context_get_double(context, "test", &value);
	CHECK_EQ(value, 0);

	fclaw_global_destroy(glob);
}

TEST_CASE("fclaw_context_get_double and fclaw_context_get_int called for same value")
{
	for(int default_value : {-100, 0, 42})
	for(double default_double : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);

		double value_double = default_double;
		fclaw_context_get_double(context, "test2", &value_double);
		CHECK_EQ(value_double, default_double);

		fclaw_context_save(context);
		
		context = fclaw_context_get(glob, "test");

		value = 0; value_double = 0;

		fclaw_context_get_int(context, "test", &value);
		fclaw_context_get_double(context, "test2", &value_double);

		CHECK_EQ(value, default_value);
		CHECK_EQ(value_double, default_double);

		fclaw_global_destroy(glob);
	}
}

TEST_CASE("fclaw_context_get_double and fclaw_context_get_int called for same key")
{
	for(int default_value : {-100, 0, 42})
	for(double default_double : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);

		double value_double = default_double;
		CHECK_SC_ABORTED(fclaw_context_get_double(context, "test", &value_double));
	}

	for(double default_double : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		double value_double = default_double;
		fclaw_context_get_double(context, "test", &value_double);

		int value;
		CHECK_SC_ABORTED(fclaw_context_get_int(context, "test", &value));
	}
}

TEST_CASE("fclaw_context pack/unpack glob int and double without save")
{
	for(int default_value : {-100, 0, 42})
	for(double default_double : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_vtable_initialize(glob);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);

		double value_double = default_double;
		fclaw_context_get_double(context, "test2", &value_double);
		CHECK_EQ(value_double, default_double);

		char buffer[fclaw_global_packsize(glob)];
		CHECK_SC_ABORTED(fclaw_global_pack(glob, buffer));
	}
}

TEST_CASE("fclaw_context pack/unpack glob int and double")
{
	for(int default_value : {-100, 0, 42})
	for(double default_double : {-100, 0, 42})
	{
		fclaw_global_t* glob = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_vtable_initialize(glob);
		fclaw_context_t *context = fclaw_context_get(glob, "test");

		int value = default_value;
		fclaw_context_get_int(context, "test", &value);
		CHECK_EQ(value, default_value);

		double value_double = default_double;
		fclaw_context_get_double(context, "test2", &value_double);
		CHECK_EQ(value_double, default_double);

		fclaw_context_save(context);

		char buffer[fclaw_global_packsize(glob)];
		fclaw_global_pack(glob, buffer);

		fclaw_global_t* glob2 = fclaw_global_new_comm(sc_MPI_COMM_SELF, 1, 0);
		fclaw_context_vtable_initialize(glob2);

		fclaw_global_unpack(buffer, glob2);

		fclaw_context_t *context2 = fclaw_context_get(glob2, "test");
		int value2 = 0;
		double value_double2 = 0;

		fclaw_context_get_int(context2, "test", &value2);
		CHECK_EQ(value2, default_value);

		fclaw_context_get_double(context2, "test2", &value_double2);
		CHECK_EQ(value_double2, default_double);

		fclaw_global_destroy(glob);
		fclaw_global_destroy(glob2);
	}
}