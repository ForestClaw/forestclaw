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

#ifndef FCLAW_SERIALIZATION_H
#define FCLAW_SERIALIZATION_H

#include <fclaw_base.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/* these are dimension-specific functions */

/**
 * @brief Pack userdata into a buffer
 * @param userdata pointer to userdata
 * @param buffer buffer to pack into
 */
typedef size_t (*fclaw_packing_pack_t)(void* userdata,
                                     char* buffer);
/**
 * @brief Unpack userdata from buffer 
 * @param buffer buffer to unpack from
 * @return newly create userdata
 */
typedef size_t (*fclaw_packing_unpack_t)(char* buffer,void**);
/**
 * @brief Get the size needed to pack userdata
 * @return the size
 */
typedef size_t (*fclaw_packing_packsize_t)(void* userdata);

/**
 * @brief destroy userdata
 */
typedef void (*fclaw_packing_destroy_t)(void* value);

/**
 * @brief vtable for packing functions
 */
typedef struct fclaw_packing_vtable
{
  fclaw_packing_pack_t pack; /**< function for packing */
  fclaw_packing_unpack_t unpack; /**< function for unpacking*/
  fclaw_packing_packsize_t size; /**< function for packing size */
  fclaw_packing_destroy_t destroy; /**< function for destroying */
} fclaw_useradata_vtable_t;

/**
 * @brief Get the size needed for packing a string
 * 
 * @return size_t number of bytes needed
 */
size_t fclaw_packsize_string(const char*);

/**
 * @brief Pack a string into a buffer
 * 
 * @param buffer the buffer
 * @return size_t bytes written
 */
size_t fclaw_pack_string(const char*, char* buffer);

/**
 * @brief Unpack a string from a buffer
 * 
 * @param buffer the buffer
 * @return size_t bytes read
 */
size_t fclaw_unpack_string(const char* buffer, char**);

/**
 * @brief Pack an integer into a buffer 
 * 
 * @param value the int 
 * @param buffer the buffer
 * @return size_t the number of bytes written
 */
size_t fclaw_pack_int(int value, char* buffer);

/**
 * @brief Unpack an integer from a buffer
 * 
 * @param buffer the buffer
 * @param value the int
 * @return size_t the number of bytes read
 */
size_t fclaw_unpack_int(const char* buffer, int* value);

/**
 * @brief Pack a size_t into buffer
 * 
 * @param value the size_t
 * @param buffer the buffer
 * @return size_t number of bytes writen
 */
size_t fclaw_pack_size_t(size_t value, char * buffer);

/**
 * @brief Unpack a size_t from buffer
 * 
 * @param buffer the buffer
 * @param value the value
 * @return size_t number of bytes read
 */
size_t fclaw_unpack_size_t(const char* buffer, size_t* value);

/**
 * @brief Pack a double into a buffer
 * 
 * @param value the double
 * @param buffer the buffer
 * @return size_t 
 */
size_t fclaw_pack_double(double value, char * buffer);

/**
 * @brief Unpack a double from a buffer
 * 
 * @param buffer the buffer
 * @param value the value  
 * @return size_t the number of bytes read
 */
size_t fclaw_unpack_double(const char * buffer, double* value);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_GLOBAL_H */
