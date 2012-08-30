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

#ifndef AMR_OPTIONS_H
#define AMR_OPTIONS_H

#include <sc_options.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

typedef struct amr_options
{
    int mx_leaf, my_leaf;
    double tfinal;
    const char *subcycling;
}
amr_options_t;

/* add options specific to forestclaw to an existing sc_options structure */
void amr_options_register (sc_options_t * opt, amr_options_t * amropt);

/* parse options and populate values in registered amr_options structure */
void amr_options_parse (sc_options_t * opt, int argc, char ** argv,
                        int log_priority);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !AMR_OPTIONS_H */
