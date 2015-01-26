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

#include <fclaw_base.h>

static void
run_program (fclaw_app_t * a)
{
    int debug_size, debug_rank;

    fclaw_global_essentialf ("So this is the beginning of the program\n");
    (void) fclaw_app_get_mpi_size_rank (a, &debug_size, &debug_rank);

    /* this is where we would do some numerics */
    fclaw_debugf ("Debug message (rank %d/%d)\n", debug_rank, debug_size);
    fclaw_infof ("Info message (individual)\n");
    fclaw_global_infof ("Info message\n");
    fclaw_global_productionf ("Production message\n");
    fclaw_global_essentialf ("Essential message\n");

    fclaw_global_essentialf ("And this is the end of the program\n");
}

int
main (int argc, char **argv)
{
#if 0
    int dummyvar;
#endif
    int first_arg;
    fclaw_exit_type_t vexit;
    fclaw_app_t *a;

    /* initialize application */
    a = fclaw_app_new (&argc, &argv, NULL);
    fclaw_app_options_register_core (a, NULL);
    /* TODO: register more options packages here */
    vexit = fclaw_app_options_parse (a, &first_arg);

#if 0
    /* THIS WILL BE DONE PER-PACKAGE IN AN INTERFACE YET TO BE DEVELOPED */
    /* register options */
    sc_options_add_int (opt, '\0', "dummy", &dummyvar, 5, "Dummy integer");
#endif

    if (!vexit)
    {
        /* parameters are clean */
        run_program (a);
    }

    /* cleanup application */
    fclaw_app_destroy (a);
    return fclaw_app_exit_type_to_status (vexit);
}
