/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW_CLAWPATCH_OUTPUT_ASCII_H
#define FCLAW_CLAWPATCH_OUTPUT_ASCII_H

#ifdef __cplusplus
extern "C"
{
#endif

struct fclaw_global;
struct fclaw_patch;
struct fclaw_domain;

/** 
 * @file
 * Routines for ascii output 
 */

/**
 * @brief Callback called on each patch
 * 
 * @param domain the domain context
 * @param patch the patch context
 * @param blockno the block number
 * @param patchno the patch number
 * @param user the user data pointer
 */
void cb_clawpatch_output_ascii (struct fclaw_domain * domain,
                                struct fclaw_patch * patch,
                                int blockno, int patchno,
                                void *user);

/**
 * @brief output ascii data
 * 
 * @param glob the global context
 * @param iframe the frame index
 */
void fclaw_clawpatch_output_ascii(struct fclaw_global* glob,int iframe);

/**
 * @brief output ascii time header
 * 
 * @param glob the global context
 * @param iframe the frame index
 */
void fclaw_clawpatch_time_header_ascii(struct fclaw_global* glob, int iframe);


#ifdef __cplusplus
}
#endif

#endif
