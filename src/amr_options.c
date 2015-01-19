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

#include <amr_options.h>


/* This is here for backwards compatibility */
amr_options_t * amr_options_new (sc_options_t * opt)
{
    amr_options_t *amropt = fclaw_options_new();
    fclaw_options_register(opt,amropt);
    sc_options_load (sc_package_id, SC_LP_ALWAYS, opt,
                     "fclaw_options.ini");
    return amropt;
}

void amr_options_parse (sc_options_t * opt,
                        int argc, char **argv, int log_priority)
{
    fclaw_options_parse_command_line(opt,argc,argv);
}

void amr_options_destroy (amr_options_t * gparms)
{
    fclaw_options_destroy_array((void*) gparms->mthbc);

    /* We shouldn't do this if the user has a static declaration of gparms */
    FCLAW_FREE(gparms);
}

void amr_postprocess_parms (amr_options_t * amropt)
{
    fclaw_options_postprocess(amropt);
}


void amr_checkparms (amr_options_t * gparms)
{
    /* Check outstyle. */
    if (gparms->outstyle == 1 && gparms->use_fixed_dt)
    {
        double dT_outer = gparms->tfinal / gparms->nout;
        double dT_inner = gparms->initial_dt;
        int nsteps = (int) floor (dT_outer / dT_inner + .5);
        if (fabs (nsteps * dT_inner - dT_outer) > 1e-8)
        {
            printf
                ("For fixed dt, initial time step size must divide tfinal/nout "
                 "exactly.\n");

            /*
             * TODO: exit should only be called from the main program.
             * In forestclaw, probably not even there:
             * We should exit gracefully by calling mpi_finalize etc.
             */
            exit (1);
        }
    }

    /* Could also do basic sanity checks on mx,my,... */
}
