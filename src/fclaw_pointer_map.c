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

#include <fclaw_base.h>
#include <fclaw_pointer_map.h>
#include <sc_keyvalue.h>

typedef struct value
{
    void* pointer;
    fclaw_pointer_map_value_destroy_t destroy;
}value_t;

fclaw_pointer_map_t* fclaw_pointer_map_new()
{
    return (fclaw_pointer_map_t*) sc_keyvalue_new();
}

static int 
destroy_iter(const char *key,
             const sc_keyvalue_entry_type_t
             type, void *entry,
             const void *u)
{
    value_t* value = *((value_t**) entry);

    if(value != NULL){
        if(value->destroy != NULL)
            value->destroy(value->pointer);

        FCLAW_FREE(value);
    }

    FCLAW_FREE((void*) key);

    return 1;
}

void fclaw_pointer_map_destroy(fclaw_pointer_map_t* map)
{
    sc_keyvalue_t *kv = (sc_keyvalue_t*) map;
    
    sc_keyvalue_foreach(kv, destroy_iter, NULL);

    sc_keyvalue_destroy(kv);
}

void fclaw_pointer_map_insert(fclaw_pointer_map_t* map, 
                              const char* key, 
                              void* pointer, 
                              fclaw_pointer_map_value_destroy_t destroy)
{
    sc_keyvalue_t *kv = (sc_keyvalue_t*) map;

    value_t* old_value = (value_t*) sc_keyvalue_get_pointer(kv, key, NULL);

    if(old_value != NULL){
        if(old_value->destroy != NULL){
            old_value->destroy(old_value->pointer);
        }
        FCLAW_FREE(old_value);
    }

    value_t* value = FCLAW_ALLOC(value_t, 1);
    value->pointer = pointer;
    value->destroy = destroy;

    char* key_copy = FCLAW_ALLOC(char, strlen(key) + 1);
    strcpy(key_copy, key);

    sc_keyvalue_set_pointer(kv, key_copy, value);
}

void* fclaw_pointer_map_get(fclaw_pointer_map_t* map, const char* key)
{
    sc_keyvalue_t *kv = (sc_keyvalue_t*) map;

    value_t* value = (value_t*) sc_keyvalue_get_pointer(kv, key, NULL);

    if(value == NULL)
        return NULL;

    return value->pointer;
}

struct userdata
{
    fclaw_pointer_map_iterate_cb_t cb;
    void* user;
};

static int 
iterator_iter(const char *key,
              const sc_keyvalue_entry_type_t
              type, void *entry,
              const void *u)
{
    struct userdata* userdata = (struct userdata*) u;
    value_t* value = *((value_t**) entry);
    userdata->cb(key, value->pointer, userdata->user);
    return 1;
}

void fclaw_pointer_map_iterate(fclaw_pointer_map_t *map, fclaw_pointer_map_iterate_cb_t cb, void *user){
    sc_keyvalue_t *kv = (sc_keyvalue_t*) map;

    struct userdata userdata;
    userdata.user = user;
    userdata.cb = cb;

    sc_keyvalue_foreach(kv, iterator_iter, &userdata);
}

static int 
size_iter(const char *key,
     const sc_keyvalue_entry_type_t
     type, void *entry,
     const void *u)
{
    *(int*)u += 1;
    return 1;
}

int fclaw_pointer_map_size(fclaw_pointer_map_t* map){
    sc_keyvalue_t *kv = (sc_keyvalue_t*) map;
    int size = 0;
    sc_keyvalue_foreach(kv, size_iter, &size);
    return size;
}