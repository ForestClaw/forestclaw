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

#ifndef AMR_OUTPUT_H
#define AMR_OUTPUT_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

typedef struct fclaw2d_output_functions fclaw2d_output_functions_t;

typedef void (*fclaw2d_patch_write_header_t)(fclaw2d_domain_t* domain, int iframe, int ngrids);

typedef void (*fclaw2d_patch_write_output_t)(fclaw2d_domain_t *domain,
                                             fclaw2d_patch_t *this_patch,
                                             int this_block_idx, int this_patch_idx,
                                             int iframe,int patch_num,int level);

/* --------------------------------------------------------------
   Some default functions - the user can write their own versions
   of these to get functionality they want.
   ----------------------------------------------------------------- */
void matlab_write_header(fclaw2d_domain_t* domain, int iframe,int ngrids);

void matlab_write_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                         int this_block_idx, int this_patch_idx,
                         int iframe, int patch_num, int level);

/* --------------------------------------------------------------
   Generic functions for writing output
   ----------------------------------------------------------------- */

fclaw2d_output_functions_t* get_output_functions(fclaw2d_domain_t *domain);

void link_output_functions(fclaw2d_domain_t* domain,
                           fclaw2d_patch_write_header_t f_patch_write_header,
                           fclaw2d_patch_write_output_t f_patch_write_output);

void initialize_output_functions(fclaw2d_output_functions_t *output_functions);

void copy_output_functions(fclaw2d_output_functions_t *old_functions,
                           fclaw2d_output_functions_t *new_functions);

/* ------------------------------------------------------------------------- */
struct fclaw2d_output_functions
{
    fclaw2d_patch_write_header_t   f_patch_write_header;
    fclaw2d_patch_write_output_t   f_patch_write_output;
};


#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
