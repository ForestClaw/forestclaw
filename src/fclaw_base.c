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

static int fclaw_package_id = -1;

int
fclaw_app_exit_type_to_status (fclaw_exit_type_t vexit)
{
    return vexit == FCLAW_EXIT_ERROR ? EXIT_FAILURE : EXIT_SUCCESS;
}

/** Each options packages lives in a structure like this. */
typedef struct fclaw_app_options
{
    char *section;              /**< NULL or used in sc_options_add_suboptions. */
    char *configfile;           /**< NULL or an .ini file's basename to read. */
    fclaw_app_options_vtable_t vt;      /**< Virtual table for option processing. */
    void *package;              /**< The package user data from options_register. */
    void *registered;           /**< Whatever is returend by options_register. */
}
fclaw_app_options_t;

/** An application container whose use is optional. */
struct fclaw_app
{
    sc_MPI_Comm mpicomm;      /**< Communicator is set to MPI_COMM_WORLD. */
    int mpisize;              /**< Size of communicator. */
    int mpirank;              /**< Rank of this process in \b mpicomm. */
    int first_arg;            /**< Location of first non-option argument after parsing. */
    int *argc;                /**< Pointer to main function's argument count. */
    char ***argv;             /**< Pointer to main function's argument list. */
    sc_options_t *opt;        /**< Central options structure. */
    void *user;               /**< Set by fclaw_app_new, not touched by forestclaw. */

    /* options packages */
    sc_array_t *opt_pkg;      /**< An array of fclaw_app_options types. */
};

int
fclaw_get_package_id (void)
{
    return fclaw_package_id;
}

static void
fclaw_logv (int category, int priority, const char *fmt, va_list ap)
{
    sc_logv ("unknown", -1, fclaw_package_id, category, priority, fmt, ap);
}

void
fclaw_logf (int category, int priority, const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    sc_logv ("unknown", -1, fclaw_package_id, category, priority, fmt, ap);
    va_end (ap);
}

void
fclaw_global_essentialf (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_GLOBAL, SC_LP_ESSENTIAL, fmt, ap);
    va_end (ap);
}

void
fclaw_global_productionf (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_GLOBAL, SC_LP_PRODUCTION, fmt, ap);
    va_end (ap);
}

void
fclaw_global_infof (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_GLOBAL, SC_LP_INFO, fmt, ap);
    va_end (ap);
}

void
fclaw_infof (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_NORMAL, SC_LP_INFO, fmt, ap);
    va_end (ap);
}

void
fclaw_debugf (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_NORMAL, SC_LP_DEBUG, fmt, ap);
    va_end (ap);
}

void
fclaw_init (sc_log_handler_t log_handler, int log_threshold)
{
    int w;

    fclaw_package_id = sc_package_register (log_handler, log_threshold,
                                            "fclaw",
                                            "Solve conservation laws");

    w = 24;
    fclaw_global_essentialf ("This is %s\n", FCLAW_PACKAGE_STRING);
    fclaw_global_productionf ("%-*s %s\n", w, "CPP", FCLAW_CPP);
    fclaw_global_productionf ("%-*s %s\n", w, "CPPFLAGS", FCLAW_CPPFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "CC", FCLAW_CC);
    fclaw_global_productionf ("%-*s %s\n", w, "CFLAGS", FCLAW_CFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "CXX", FCLAW_CXX);
    fclaw_global_productionf ("%-*s %s\n", w, "CXXFLAGS", FCLAW_CXXFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "LDFLAGS", FCLAW_LDFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "LIBS", FCLAW_LIBS);
}

fclaw_app_t *
fclaw_app_new (int *argc, char ***argv, void *user)
{
#ifdef FCLAW_ENABLE_DEBUG
    const int LP_lib = SC_LP_INFO;
    const int LP_fclaw = SC_LP_DEBUG;
#else
    const int LP_lib = SC_LP_ESSENTIAL;
    const int LP_fclaw = SC_LP_PRODUCTION;
#endif
    int mpiret;
    MPI_Comm mpicomm;
    fclaw_app_t *a;

    mpiret = sc_MPI_Init (argc, argv);
    SC_CHECK_MPI (mpiret);
    mpicomm = sc_MPI_COMM_WORLD;

    sc_init (mpicomm, 1, 1, NULL, LP_lib);
    p4est_init (NULL, LP_lib);
    fclaw_init (NULL, LP_fclaw);

    a = FCLAW_ALLOC (fclaw_app_t, 1);
    a->mpicomm = mpicomm;
    mpiret = sc_MPI_Comm_size (a->mpicomm, &a->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (a->mpicomm, &a->mpirank);
    SC_CHECK_MPI (mpiret);

    srand (a->mpirank);
    a->first_arg = -1;
    a->argc = argc;
    a->argv = argv;
    a->user = user;
    a->opt = sc_options_new ((*argv)[0]);
    a->opt_pkg = sc_array_new (sizeof (fclaw_app_options_t));

    return a;
}

void
fclaw_app_destroy (fclaw_app_t * a)
{
    int mpiret;
    size_t zz;
    fclaw_app_options_t *ao;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (a->opt_pkg != NULL);
    FCLAW_ASSERT (a->opt != NULL);

    /* let the options packages clean up their memory */
    for (zz = 0; zz < a->opt_pkg->elem_count; ++zz)
    {
        ao = (fclaw_app_options_t *) sc_array_index (a->opt_pkg, zz);
        FCLAW_ASSERT (ao != NULL);
        if (ao->vt.options_destroy != NULL)
        {
            ao->vt.options_destroy (a, ao->package, ao->registered);
        }
    }
    sc_array_destroy (a->opt_pkg);

    /* free central structures */
    sc_options_destroy (a->opt);
    FCLAW_FREE (a);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
}

void
fclaw_app_options_register (fclaw_app_t * a,
                            const char *section, const char *configfile,
                            const fclaw_app_options_vtable_t * vt,
                            void *package)
{
    sc_options_t *popt;
    fclaw_app_options_t *ao;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (vt != NULL && vt->options_register != NULL);

    ao = (fclaw_app_options_t *) sc_array_push (a->opt_pkg);
    ao->section = section == NULL ? NULL : FCLAW_STRDUP (section);
    ao->configfile = configfile == NULL ? NULL : FCLAW_STRDUP (section);
    ao->vt = *vt;
    ao->package = package;

    popt = section == NULL ? a->opt : sc_options_new (section);
    ao->registered = vt->options_register (a, package, popt);
    if (section != NULL)
    {
        FCLAW_ASSERT (popt != a->opt);
        sc_options_add_suboptions (a->opt, popt, section);
        sc_options_destroy (popt);
    }
}

fclaw_exit_type_t
fclaw_app_options_parse (fclaw_app_t * a, int * first_arg)
{
  return FCLAW_NOEXIT;
}

MPI_Comm
fclaw_app_get_mpi_size_rank (fclaw_app_t * a, int *mpisize, int *mpirank)
{
    FCLAW_ASSERT (a != NULL);

    if (mpisize != NULL)
    {
        *mpisize = a->mpisize;
    }
    if (mpirank != NULL)
    {
        *mpirank = a->mpirank;
    }
    return a->mpicomm;
}

void *
fclaw_app_get_user (fclaw_app_t * a)
{
    FCLAW_ASSERT (a != NULL);

    return a->user;
}

sc_options_t *
fclaw_app_get_options (fclaw_app_t * a)
{
    FCLAW_ASSERT (a != NULL);

    return a->opt;
}
