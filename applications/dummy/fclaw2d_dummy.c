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

static int fclaw_package_id;

static void
run_program (fclaw_app_t * a)
{
    fclaw_global_essentialf ("So this is the beginning of the program\n");

    /* this is where we would do some numerics */
    fclaw_debugf ("Debug message (individual)\n");
    fclaw_infof ("Info message (individual)\n");
    fclaw_global_infof ("Info message\n");
    fclaw_global_productionf ("Production message\n");
    fclaw_global_essentialf ("Essential message\n");

    fclaw_global_essentialf ("And this is the end of the program\n");
}

int
main (int argc, char **argv)
{
    int dummyvar;
    int verbosity;
    int first_arg;
    sc_keyvalue_t *kv_verbosity;
    fclaw_app_t sapp, *a = &sapp;

    /* initialize application */
    fclaw_app_init (a, &argc, &argv, NULL);
    fclaw_package_id = fclaw_get_package_id ();

    /* THIS WILL BE DONE PER-PACKAGE IN AN INTERFACE YET TO BE DEVELOPED */
    /* register options */
    sc_options_add_int (a->opt, '\0', "dummy", &dummyvar, 5, "Dummy integer");
    kv_verbosity = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv_verbosity, "default", FCLAW_VERBOSITY_DEFAULT);
    sc_keyvalue_set_int (kv_verbosity, "debug", FCLAW_VERBOSITY_DEBUG);
    sc_keyvalue_set_int (kv_verbosity, "info", FCLAW_VERBOSITY_INFO);
    sc_keyvalue_set_int (kv_verbosity, "production",
                         FCLAW_VERBOSITY_PRODUCTION);
    sc_keyvalue_set_int (kv_verbosity, "essential",
                         FCLAW_VERBOSITY_ESSENTIAL);
    sc_keyvalue_set_int (kv_verbosity, "silent", FCLAW_VERBOSITY_SILENT);
    sc_options_add_keyvalue (a->opt, 'V', "fclaw-verbosity", &verbosity,
                             "default", kv_verbosity, "Set verbosity level");

    /* THIS BLOCK WILL BE ENCAPSULATED -- NOT FOR APPLICATION EXAMPLES */
    /* we're still missing to load/save .ini files */
    /* parse command line options */
    first_arg = sc_options_parse (sc_package_id, FCLAW_VERBOSITY_ESSENTIAL,
                                  a->opt, argc, argv);
    if (first_arg < 0)
    {
        fclaw_global_essentialf ("Option parsing failed\n");
        sc_options_print_usage (fclaw_package_id, FCLAW_VERBOSITY_INFO,
                                a->opt, NULL);
        fclaw_global_essentialf ("Terminating program\n");
    }
    else
    {
        fclaw_global_infof ("Option parsing successful\n");
        sc_options_print_summary (fclaw_package_id,
                                  FCLAW_VERBOSITY_PRODUCTION, a->opt);

        /* set verbosity levels */
        sc_package_set_verbosity (sc_package_id, FCLAW_VERBOSITY_ESSENTIAL);
        sc_package_set_verbosity (p4est_package_id,
                                  FCLAW_VERBOSITY_ESSENTIAL);
        sc_package_set_verbosity (fclaw_package_id, verbosity);

        /* go to work */
        run_program (a);
    }

    /* cleanup application */
    fclaw_app_reset (a);

    return 0;
}
