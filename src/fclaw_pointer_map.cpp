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
#include <map>
#include <string>

namespace
{
    struct value
    {
        void* pointer = nullptr;
        fclaw_pointer_map_value_destroy_t destroy = nullptr;
    };
}

struct fclaw_pointer_map
{
    std::map<std::string, value> map;
};

fclaw_pointer_map_t* fclaw_pointer_map_new()
{
    return new fclaw_pointer_map();
}

void fclaw_pointer_map_destroy(fclaw_pointer_map_t* map)
{
    std::map<std::string, value>::iterator it = map->map.begin();
    while(it != map->map.end()){
        if(it->second.destroy != nullptr)
            it->second.destroy(it->second.pointer);
        ++it;
    }
    delete map;
}

void fclaw_pointer_map_insert(fclaw_pointer_map_t* map, 
                              const char* key, 
                              void* pointer, 
                              fclaw_pointer_map_value_destroy_t destroy)
{
    value& v = map->map[key];
    if(v.destroy){
        v.destroy(v.pointer);
    }
    v.pointer = pointer;
    v.destroy = destroy;
}

void* fclaw_pointer_map_get(fclaw_pointer_map_t* map, const char* key)
{
    return map->map[key].pointer;
}

