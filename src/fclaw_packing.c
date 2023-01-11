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
#include <fclaw_base.h>


size_t fclaw_packsize_string(const char * string){
  return sizeof(size_t) + (string == NULL ? 0 : strlen(string)+1);
}

size_t fclaw_pack_string(const char* string, char* buffer){
  char * buffer_start = buffer;
  if(string == NULL){
    buffer += fclaw_pack_size_t(buffer, 0);
  }else{
    size_t length = strlen(string)+1;
    buffer += fclaw_pack_size_t(buffer, length);
    memcpy(buffer,string,length);
    buffer += length;
  }
  return buffer - buffer_start;
}

size_t fclaw_unpack_string(char * buffer, char** string){
  char * buffer_start = buffer;
  size_t length;
  buffer += fclaw_unpack_size_t(buffer, &length);
  if(length == 0){
    *string = NULL;
  }else{
    *string = FCLAW_ALLOC(char,length);
    memcpy(*string,buffer,length);
    buffer += length;
  }
  return buffer - buffer_start;
}

size_t fclaw_pack_int(int value, char * buffer){
  *((int *) buffer) = value;
  return sizeof(int);
}

size_t fclaw_unpack_int(char * buffer, int* value){
  *value = *((int *) buffer);
  return sizeof(int);
}

size_t fclaw_pack_size_t(char * buffer, size_t value){
  *((size_t *) buffer) = value;
  return sizeof(size_t);
}

size_t fclaw_unpack_size_t(char * buffer, size_t* value){
  *value = *((size_t *) buffer);
  return sizeof(size_t);
}

size_t fclaw_pack_double(double value, char* buffer){
  *((double *) buffer) = value;
  return sizeof(double);
}

size_t fclaw_unpack_double(const char * buffer, double* value){
  *value = *((double *) buffer);
  return sizeof(double);
}