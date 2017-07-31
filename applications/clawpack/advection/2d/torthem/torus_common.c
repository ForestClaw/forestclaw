/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
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

#include "torus_common.h"

fclaw_exit_type_t
torus_options_postprocess (user_options_t * user)
{
    if (user->example == 3)
    {
        fclaw_options_convert_double_array (user->latitude_string,
                                            &user->latitude, 2);
        fclaw_options_convert_double_array (user->longitude_string,
                                            &user->longitude, 2);
    }
    return FCLAW_NOEXIT;
}

fclaw_exit_type_t
torus_options_check (user_options_t * user)
{
    if (user->example < 0 || user->example > 7)
    {
        fclaw_global_essentialf
            ("Option --user:example must be between 0 and 7\n");
        return FCLAW_EXIT_QUIET;
    }
    return FCLAW_NOEXIT;
}

void
torus_options_destroy (user_options_t * user)
{
    fclaw_options_destroy_array (user->latitude);
    fclaw_options_destroy_array (user->longitude);
}
