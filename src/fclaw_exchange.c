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

#include <sc.h>
#include <fclaw_global.h>
#include <fclaw_domain.h>
#include <fclaw_options.h>
#include <fclaw_patch.h>
#include <fclaw_exchange.h>

// 2d/3d functions hidden from users
void fclaw2d_exchange_setup(struct fclaw_global* glob,
                            fclaw_timer_names_t running);

void fclaw2d_exchange_delete(struct fclaw_global* glob);

void fclaw2d_exchange_ghost_patches_begin(struct fclaw_global* glob,
                                          int minlevel,
                                          int maxlevel,
                                          int time_interp,
                                          fclaw_timer_names_t running);

void fclaw2d_exchange_ghost_patches_end(struct fclaw_global* glob,
                                        int minlevel,
                                        int maxlevel,
                                        int time_interp,
                                        fclaw_timer_names_t running);
void fclaw3d_exchange_setup(struct fclaw_global* glob,
                            fclaw_timer_names_t running);

void fclaw3d_exchange_delete(struct fclaw_global* glob);

void fclaw3d_exchange_ghost_patches_begin(struct fclaw_global* glob,
                                          int minlevel,
                                          int maxlevel,
                                          int time_interp,
                                          fclaw_timer_names_t running);

void fclaw3d_exchange_ghost_patches_end(struct fclaw_global* glob,
                                        int minlevel,
                                        int maxlevel,
                                        int time_interp,
                                        fclaw_timer_names_t running);
/* --------------------------------------------------------------------------
   Public interface
   -------------------------------------------------------------------------- */
/* This is called whenever a new domain is created (initialize, regrid) */
void fclaw_exchange_setup(fclaw_global_t* glob,
                            fclaw_timer_names_t running)
{
    if(glob->domain->dim == 2)
    {
        fclaw2d_exchange_setup(glob,running);
    }
    else if (glob->domain->dim == 3)
    {
        fclaw3d_exchange_setup(glob,running);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }

}

void fclaw_exchange_delete(fclaw_global_t* glob)
{
    if(glob->domain->dim == 2)
    {
        fclaw2d_exchange_delete(glob);
    }
    else if (glob->domain->dim == 3)
    {
        fclaw3d_exchange_delete(glob);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}


/* ----------------------------------------------------------------
   Public interface
   -------------------------------------------------------------- */

/* This is called whenever all time levels are time synchronized. */
void fclaw_exchange_ghost_patches_begin(fclaw_global_t* glob,
                                          int minlevel,
                                          int maxlevel,
                                          int time_interp,
                                          fclaw_timer_names_t running)
{
    if(glob->domain->dim == 2)
    {
        fclaw2d_exchange_ghost_patches_begin(glob,minlevel,maxlevel,
                                             time_interp,running);
    }
    else if (glob->domain->dim == 3)
    {
        fclaw3d_exchange_ghost_patches_begin(glob,minlevel,maxlevel,
                                             time_interp,running);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}

/* This is called whenever all time levels are time synchronized. */
void fclaw_exchange_ghost_patches_end(fclaw_global_t* glob,
                                        int minlevel,
                                        int maxlevel,
                                        int time_interp,
                                        fclaw_timer_names_t running)
{
    if(glob->domain->dim == 2)
    {
        fclaw2d_exchange_ghost_patches_end(glob,minlevel,maxlevel,
                                           time_interp,running);
    }
    else if (glob->domain->dim == 3)
    {
        fclaw3d_exchange_ghost_patches_end(glob,minlevel,maxlevel,
                                           time_interp,running);
    }
    else
    {
        SC_ABORT_NOT_REACHED();
    }
}
