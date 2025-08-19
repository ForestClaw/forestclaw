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

#include <fclaw_context.h>
#include <fclaw_packing.h>

#define PACKING_VTABLE_NAME "fclaw_context_t"

typedef enum
{
    FCLAW_CONTEXT_INT,
    FCLAW_CONTEXT_DOUBLE,
    FCLAW_CONTEXT_POINTER
} type_t;

typedef struct value
{
    type_t type;
    union
    {
        int i;
        double d;
    } value;

    int initializing;

    void *pointer;

} value_t;

static void
value_destroy(void *data)
{
    value_t *value = (value_t *)data;
    FCLAW_FREE(value);
}

struct fclaw_context
{
    int saved;
    fclaw_pointer_map_t *values;
};

static fclaw_context_t *
context_new()
{
    fclaw_context_t *context = FCLAW_ALLOC(fclaw_context_t, 1);
    context->saved = 0;
    context->values = fclaw_pointer_map_new();
    return context;
}

static void 
context_destroy(void *data)
{
    fclaw_context_t *context = (fclaw_context_t *)data;
    fclaw_pointer_map_destroy(context->values);
    FCLAW_FREE(context);
}

static void
reset_values(const char *key, void *data, void *user)
{
    value_t *value = (value_t *)data;
    value->pointer = NULL;
    value->initializing = 1;
}

fclaw_context_t* fclaw_context_get(fclaw_global_t *glob, const char *name)
{
    fclaw_context_t *context = 
        (fclaw_context_t*) fclaw_global_get_attribute(glob, name);

    if(context == NULL)
    {
        context = context_new();
        fclaw_global_attribute_store(glob, name, context, PACKING_VTABLE_NAME, context_destroy);
    }
    else
    {
        if(!context->saved)
        {
            fclaw_abortf("fclaw_context_get: Context needs to be saved before it can be retrieved again\n");
        }
        context->saved = 0;
        fclaw_pointer_map_iterate(context->values, reset_values, NULL);
    }

    return context;
}


void fclaw_context_get_int(fclaw_context_t *context, 
                           const char *name,
                           int *value)
{
    value_t *v = (value_t*) fclaw_pointer_map_get(context->values, name);
    if(v != NULL)
    {
        if(v->type != FCLAW_CONTEXT_INT)
        {
            fclaw_abortf("fclaw_context_get_int: Value %s is not an int\n", name);
        }
        if(v->initializing)
        {
            *value = v->value.i;
            v->initializing = 0;
        }
    }
    else
    {
        v = FCLAW_ALLOC(value_t, 1);
        v->type = FCLAW_CONTEXT_INT;
        v->initializing = 0;
        fclaw_pointer_map_insert(context->values, name, v, value_destroy);
    }
    
    v->pointer = value;
}

void fclaw_context_get_double(fclaw_context_t *context, 
                              const char *name, 
                              double *value)
{
    value_t *v = (value_t*) fclaw_pointer_map_get(context->values, name);
    if(v != NULL)
    {
        if(v->type != FCLAW_CONTEXT_DOUBLE)
        {
            fclaw_abortf("fclaw_context_get_double: Value %s is not a double\n", name);
        }
        if(v->initializing)
        {
            *value = v->value.d;
            v->initializing = 0;
        }
    }
    else
    {
        v = FCLAW_ALLOC(value_t, 1);
        v->type = FCLAW_CONTEXT_DOUBLE;
        v->initializing = 0;
        fclaw_pointer_map_insert(context->values, name, v, value_destroy);
    }

    v->pointer = value;
}

static
void save_value(const char *key, void *data, void *user)
{
    value_t *value = (value_t *) data;
    if(value->pointer == NULL)
    {
        /* Pointer is not set, nothing to save, old existing value will be kept */
        return;
    }
    if(value->type == FCLAW_CONTEXT_INT)
    {
        value->value.i = *(int *)value->pointer;
    }
    else if(value->type == FCLAW_CONTEXT_DOUBLE)
    {
        value->value.d = *(double *)value->pointer;
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

void fclaw_context_save(fclaw_context_t *context)
{
    fclaw_pointer_map_iterate(context->values, save_value, NULL);
    context->saved = 1;
}


static
void pack_value(const char *key, void *data, void *user)
{
    char **buffer = (char **)user;
    value_t *value = (value_t *)data;
    *buffer += fclaw_pack_string(key, *buffer);
    *buffer += fclaw_pack_int(value->type, *buffer);
    if(value->type == FCLAW_CONTEXT_INT)
    {
        *buffer += fclaw_pack_int(value->value.i, *buffer);
    }
    else if(value->type == FCLAW_CONTEXT_DOUBLE)
    {
        *buffer += fclaw_pack_double(value->value.d, *buffer);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

static
size_t pack_context(fclaw_global_t *glob, void *data, char *buffer)
{
    char* buffer_start = buffer;
    fclaw_context_t *context = (fclaw_context_t *)data;

    if(!context->saved)
    {
        fclaw_abortf("fclaw_context: Context not saved, cannot pack\n");
    }

    buffer += fclaw_pack_int(fclaw_pointer_map_size(context->values), buffer);
    fclaw_pointer_map_iterate(context->values, pack_value, &buffer);
    return buffer - buffer_start;
}

static
void value_packsize(const char *key, void *data, void *user)
{
    size_t *size = (size_t *)user;
    value_t *value = (value_t *)data;
    *size += fclaw_packsize_string(key);
    *size += sizeof(int);
    if(value->type == FCLAW_CONTEXT_INT)
    {
        *size += sizeof(int);
    }
    else if(value->type == FCLAW_CONTEXT_DOUBLE)
    {
        *size += sizeof(double);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

static
size_t context_packsize(fclaw_global_t *glob, void *data)
{
    fclaw_context_t *context = (fclaw_context_t *)data;
    size_t size = sizeof(int);
    fclaw_pointer_map_iterate(context->values, value_packsize, &size);
    return size;
}

static
size_t context_unpack(fclaw_global_t *glob, char *buffer, void *data)
{
    char* buffer_start = buffer;
    fclaw_context_t *context = (fclaw_context_t *)data;
    context->saved = 1;
    int size;
    buffer += fclaw_unpack_int(buffer, &size);
    int i;
    for(i = 0; i < size; ++i)
    {
        char *key;
        buffer += fclaw_unpack_string(buffer, &key);
        value_t *value = FCLAW_ALLOC(value_t, 1);
        buffer += fclaw_unpack_int(buffer,(int*) &value->type);
        if(value->type == FCLAW_CONTEXT_INT)
        {
            buffer += fclaw_unpack_int(buffer, &value->value.i);
        }
        else if(value->type == FCLAW_CONTEXT_DOUBLE)
        {
            buffer += fclaw_unpack_double(buffer, &value->value.d);
        }
        else
        {
            SC_ABORT_NOT_REACHED();
        }
        value->initializing = 1;
        value->pointer = NULL;
        fclaw_pointer_map_insert(context->values, key, value, value_destroy);
        FCLAW_FREE(key);
    }
    return buffer - buffer_start;
}

static void* context_new_void(fclaw_global_t *glob)
{
    return context_new();
}

fclaw_packing_vtable_t fclaw_context_vtable = {
    pack_context,
    context_unpack,
    context_packsize,
    context_new_void,
    context_destroy
};

void fclaw_context_vtable_initialize(fclaw_global_t *glob)
{
    fclaw_global_vtable_store(glob, PACKING_VTABLE_NAME, &fclaw_context_vtable, NULL);
}