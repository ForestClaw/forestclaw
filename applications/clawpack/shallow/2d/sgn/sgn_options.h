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

#ifndef SGN_USER_H
#define SGN_USER_H

#include <fclaw2d_include_all.h>
#include <fc2d_thunderegg_fort.h>  /* For virtualized functions */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

#if 0
    /* Fix syntax highlighting */
#endif    

typedef struct sgn_options
{
    int example;
    double g;
    double a;
    double b;
    double h0;

    double breaking;
    double alpha;
    double dry_tolerance;
    double sea_level;

    int claw_version;

    int is_registered;
} sgn_options_t;


sgn_options_t* sgn_options_register (struct fclaw_app * app,
                                     const char *configfile);

void sgn_options_store (struct fclaw2d_global* glob, 
                        struct sgn_options* sgn_opt);

sgn_options_t* sgn_get_options(struct fclaw2d_global* glob);



#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
