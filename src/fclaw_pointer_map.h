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

/**
 * @file
 * 
 * @brief c interface to c++ std::map<std::string,void*>
 * 
 */
#ifndef FCLAW_POINTER_MAP_H
#define FCLAW_POINTER_MAP_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/**
 * @brief Pointer map structure
 */
typedef struct fclaw_pointer_map fclaw_pointer_map_t;

/**
 * @brief function that destroys stored value
 */
typedef void (*fclaw_pointer_map_value_destroy_t)(void* value);

/**
 * @brief create a pointer map
 * 
 * @return fclaw_pointer_map_t* the new pointer map
 */
fclaw_pointer_map_t* fclaw_pointer_map_new();

/**
 * @brief destroy a pointer map
 * 
 * @param map the map to destroy
 */
void fclaw_pointer_map_destroy(fclaw_pointer_map_t* map);

/**
 * @brief insert a value into the map
 * 
 * This will destroy the previous value for the key
 * 
 * @param map the map
 * @param key the key
 * @param value the value
 * @param destroy callback function to destroy the value
 */
void fclaw_pointer_map_insert(fclaw_pointer_map_t* map, const char* key, void* value, fclaw_pointer_map_value_destroy_t destroy);

/**
 * @brief get a value from the map
 * 
 * @param map the map
 * @param key the key
 * @return void* the value, or nullptr if not found
 */
void* fclaw_pointer_map_get(fclaw_pointer_map_t* map, const char* key);

/**
 * @brief callback for pointer map iterator
 * 
 * @param key the key
 * @param value the value associated with the key
 * @param user the user pointer passed into iterator call
 */
typedef void (*fclaw_pointer_map_iterate_cb_t)(const char* key, void* pointer, void* user);

/**
 * @brief Iterate over all key-value pairs in the map
 * 
 * @param map the map to iterate over
 * @param cb the callback called on eacy key-value pair
 * @param user user pointer passed to the callback
 */
void fclaw_pointer_map_iterate(fclaw_pointer_map_t* map,fclaw_pointer_map_iterate_cb_t cb, void* user);

/**
 * @brief Get the number of key-value pairs in the map
 * 
 * @param map the map
 * @return int the number of key-value pairs
 */
int fclaw_pointer_map_size(fclaw_pointer_map_t* map);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !FCLAW_MATH_H */
