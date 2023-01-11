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

#include <fclaw_packing.h>
#include <test.hpp>

TEST_CASE("fclaw_pack_int pack and unpack")
{
	char base_buffer[2*sizeof(int)];
	for(char* buffer : {base_buffer,base_buffer+1,base_buffer+2})
	for(int value : {-1,0,2})
	{

		int unpacked_value;

		CHECK_EQ(fclaw_pack_int(value, buffer), sizeof(int));
		CHECK_EQ(fclaw_unpack_int(buffer, &unpacked_value), sizeof(int));

		CHECK_EQ(value, unpacked_value);
	}
}

TEST_CASE("fclaw_pack_size_t pack and unpack")
{
	char base_buffer[2*sizeof(size_t)];
	for(char* buffer : {base_buffer,base_buffer+1,base_buffer+2})
	for(size_t value : {1,0,2})
	{

		size_t unpacked_value;

		CHECK_EQ(fclaw_pack_size_t(buffer,value), sizeof(size_t));
		CHECK_EQ(fclaw_unpack_size_t(buffer,&unpacked_value), sizeof(size_t));

		CHECK_EQ(value, unpacked_value);
	}
}

TEST_CASE("fclaw_pack_double pack and unpack")
{
	char base_buffer[2*sizeof(double)];
	for(char* buffer : {base_buffer,base_buffer+1,base_buffer+2})
	for(double value : {-1.3,0.2,2.9})
	{

		double unpacked_value;

		CHECK_EQ(fclaw_pack_double(value, buffer), sizeof(double));
		CHECK_EQ(fclaw_unpack_double(buffer,&unpacked_value), sizeof(double));

		CHECK_EQ(value, unpacked_value);
	}
}
TEST_CASE("fclaw_pack_string pack and unpack")
{
	char base_buffer[20];
	for(char* buffer : {base_buffer,base_buffer+1,base_buffer+2})
	for(const char* value : {"", "a", "abc"})
	{
		char* unpacked_value;

		
		CHECK_EQ(fclaw_packsize_string(value), sizeof(size_t)+1+strlen(value));
		CHECK_EQ(fclaw_pack_string(value, buffer), sizeof(size_t)+1+strlen(value));
		CHECK_EQ(fclaw_unpack_string(buffer,&unpacked_value), sizeof(size_t)+1+strlen(value));

		CHECK_EQ(strcmp(value, unpacked_value), 0);
		FCLAW_FREE(unpacked_value);
	}
}
TEST_CASE("fclaw_pack_string pack and unpack NULL")
{
	char base_buffer[20];
	for(char* buffer : {base_buffer,base_buffer+1,base_buffer+2})
	{
		char* unpacked_value;

		CHECK_EQ(fclaw_packsize_string(NULL), sizeof(size_t));
		CHECK_EQ(fclaw_pack_string(NULL, buffer), sizeof(size_t));
		CHECK_EQ(fclaw_unpack_string(buffer,&unpacked_value), sizeof(size_t));

		CHECK_EQ(unpacked_value, nullptr);
	}
}