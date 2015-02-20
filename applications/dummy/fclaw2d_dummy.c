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

/** \file fclaw2d_dummy.c
 * This file is a dummy application to showcase option and configuration handling.
 */

#include <fclaw_base.h>
#include "dummy_blackbox.h"

/** This object stores configuration variables specific to this application. */
typedef struct fclaw_dummy_options
{
    int options_int;            /**< Some integer. */
    double options_double;      /**< Some double. */
    const char *dummy_string;   /**< Some string. */
    char *dummy_storage;        /**< We demonstrate memory in an options package. */
}
fclaw_dummy_options_t;

/** Callback for registering options.
 * In this application we are using two options packages with different sections.
 * \param [in,out] a            Opaque application object.
 * \param [in,out] package      Provided to \ref fclaw_app_options_register in \ref main.
 * \param [in,out] options      We are free to add our own options to this structure.
 * \return                      What we return is available to other callbacks for
 *                              this options package.  Here we do not use the feature.
 */
static void *
fclaw_options_register (fclaw_app_t * a, void *package,
                        sc_options_t * options)
{
    fclaw_dummy_options_t *dumo = (fclaw_dummy_options_t *) package;

    sc_options_add_int (options, 'i', "integer", &dumo->options_int,
                        5, "Option integer");
    sc_options_add_double (options, 'd', "double", &dumo->options_double,
                           3., "Option double");

    return NULL;
}

/** Another callback for registering options.
 * In this application we are using two options packages with different sections.
 * This one is registered with its own section "Dummy" in a configuration file.
 * \param [in,out] a            Opaque application object.
 * \param [in,out] package      Provided to \ref fclaw_app_options_register in \ref main.
 * \param [in,out] options      We are free to add our own options to this structure.
 * \return                      We return the \b dummy_storage member of
 *                              the \ref fclaw_dummy_options structure, just for fun.
 */
static void *
fclaw_dummy_register (fclaw_app_t * a, void *package, sc_options_t * options)
{
    fclaw_dummy_options_t *dumo = (fclaw_dummy_options_t *) package;

    sc_options_add_string (options, 's', "string", &dumo->dummy_string,
                           "Dummy string value", "Dummy string");

    dumo->dummy_storage = FCLAW_STRDUP ("wonderful");

    /* We're returning this just for verification purposes */
    return dumo->dummy_storage;
}

/** Callback for cleaning up memory that we allocated in options registration.
 * This one pertains to the one registered with its own section "Dummy".
 * \param [in,out] a            Opaque application object.
 * \param [in,out] package      Provided to \ref fclaw_app_options_register in \ref main.
 * \param [in,out] registered   This is the return value to \ref fclaw_dummy_register.
 */
static void
fclaw_dummy_destroy (fclaw_app_t * a, void *package, void *registered)
{
    fclaw_dummy_options_t *dumo = (fclaw_dummy_options_t *) package;

    FCLAW_ASSERT (registered != NULL &&
                  !strcmp ((const char *) registered, "wonderful"));

    FCLAW_FREE (dumo->dummy_storage);
}

/** Virtual table for the first options package in this program. */
static const fclaw_app_options_vtable_t options_vt = {
    fclaw_options_register, NULL, NULL, NULL
};

/** Virtual table for the second options package in this program. */
static const fclaw_app_options_vtable_t dummy_vt = {
    fclaw_dummy_register, NULL, NULL,
    fclaw_dummy_destroy
};

/** Whatever the program would really do in an application.
 * \param [in,out] a            Opaque application object.  We might access it to
 *                              grab its MPI communicator or user data.
 * \param [in] dumo             The configuration values we keep for this program.
 * \param [in,out] bbox         The blackbox object.
 */
static void
run_program (fclaw_app_t * a, fclaw_dummy_options_t * dumo)
{
    int debug_size, debug_rank;
    dummy_blackbox_t * bbox;

    /* try out the application-attribute feature */
    bbox = (dummy_blackbox_t *) fclaw_app_get_attribute (a, "blackbox", NULL);
    FCLAW_ASSERT (bbox != NULL);

    fclaw_global_essentialf ("So this is the beginning of the program\n");
    (void) fclaw_app_get_mpi_size_rank (a, &debug_size, &debug_rank);

    /* this is where we would do some numerics */
    fclaw_debugf ("Debug message (rank %d/%d)\n", debug_rank, debug_size);
    fclaw_infof ("Info message (individual)\n");
    fclaw_global_infof ("Info message double is %g\n", dumo->options_double);
    fclaw_global_productionf ("Production message integer is %d\n",
                              dumo->options_int);
    fclaw_global_essentialf ("Essential message string is \"%s\"\n",
                             dumo->dummy_string);
    fclaw_global_errorf ("We demonstrate an error message\n");

    /* demonstrate the blackbox package */
    fclaw_global_productionf ("We have multiplied %d into %d\n",
                              dumo->options_int,
                              dummy_blackbox_multiply
                              (bbox, dumo->options_int));

    /* this is the end of the official work of the program */
    fclaw_global_essentialf ("And this is the end of the %s program\n",
                             dumo->dummy_storage);
}

/** Main function for the dummy application. */
int
main (int argc, char **argv)
{
    int first_arg;
    fclaw_exit_type_t vexit;
    fclaw_app_t *a;
    fclaw_dummy_options_t dummy_options, *dumo = &dummy_options;
    dummy_blackbox_t *bbox;

    /* initialize application */
    a = fclaw_app_new (&argc, &argv, NULL);

    /* this application registers the core package and two of its own. */
    fclaw_app_options_register_core (a, NULL);
    fclaw_app_options_register (a, NULL, NULL, &options_vt, dumo);
    fclaw_app_options_register (a, "Dummy", NULL, &dummy_vt, dumo);

    /* initialize a package that registers its own options package */
    bbox = dummy_blackbox_new_register (a, 4);
    fclaw_app_set_attribute (a, "blackbox", bbox);

    /* process command line options and configuration files */
    vexit = fclaw_app_options_parse (a, &first_arg, "dummy_config.ini");

    if (!vexit)
    {
        /* parameters are clean */
        run_program (a, dumo);
    }

    /* cleanup packages and application */
    fclaw_app_destroy (a);

    return fclaw_app_exit_type_to_status (vexit);
}
