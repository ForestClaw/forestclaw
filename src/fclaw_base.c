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

void
fclaw_app_init (fclaw_app_t * a, int *argc, char ***argv, void *user)
{
#ifdef FCLAW_ENABLE_DEBUG
    const int LP_lib = SC_LP_INFO;
    const int LP_fclaw = SC_LP_DEBUG;
#else
    const int LP_lib = SC_LP_ESSENTIAL;
    const int LP_fclaw = SC_LP_PRODUCTION;
#endif
    int mpiret;

    mpiret = sc_MPI_Init (argc, argv);
    SC_CHECK_MPI (mpiret);

    a->mpicomm = sc_MPI_COMM_WORLD;
    mpiret = sc_MPI_Comm_size (a->mpicomm, &a->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (a->mpicomm, &a->mpirank);
    SC_CHECK_MPI (mpiret);

    srand (a->mpirank);
    a->argc = argc;
    a->argv = argv;
    a->user = user;

    sc_init (a->mpicomm, 1, 1, NULL, LP_lib);
    p4est_init (NULL, LP_lib);
    fclaw_init (NULL, LP_fclaw);

    a->opt = sc_options_new ((*argv)[0]);
}

void
fclaw_app_reset (fclaw_app_t * a)
{
    int mpiret;

    sc_options_destroy (a->opt);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
}
