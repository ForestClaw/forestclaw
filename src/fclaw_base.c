/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_pointer_map.h>
#include <fclaw_base.h>
#include <fclaw_mpi.h>

static const char *fclaw_configdir = ".forestclaw";
static const char *fclaw_env_configdir = "FCLAW_INI_DIR";
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
    int core_registered;    /**< True if core options (-h,etc) have been registered */

    /* paths and configuration files */
    const char *configdir;    /**< Defaults to fclaw_configdir under $HOME, may
                                   be changed with \ref fclaw_app_set_configdir. */
    const char *env_configdir;          /**< Name of environment variable for a
                                             directory to find configuration files.
                                             Defaults to fclaw_env_configdir. */

    /* options packages */
    sc_array_t *opt_pkg;      /**< An array of fclaw_app_options types. */

    /* attributes */
    sc_keyvalue_t *attributes;  /**< The storage for application attributes. */
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
fclaw_global_errorf (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_GLOBAL, SC_LP_ERROR, fmt, ap);
    va_end (ap);
}

void
fclaw_errorf (const char *fmt, ...)
{
    va_list ap;

    va_start (ap, fmt);
    fclaw_logv (SC_LC_NORMAL, SC_LP_ERROR, fmt, ap);
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

static int logging_rank = 0;
static const char* logging_prefix = NULL;

void 
fclaw_set_logging_prefix(const char* new_name)
{
    logging_prefix=new_name;
}

static void
log_handler (const char *name, FILE * log_stream, const char *filename, int lineno,
                int package, int category, int priority, const char *msg)
{
    int                 wi = 0;
    int                 lindent = 0;

    wi = (category == SC_LC_NORMAL);

    if(logging_prefix != NULL){
        fprintf(log_stream, "[%s]",logging_prefix);
    }
    fputc ('[', log_stream);
    fprintf (log_stream, "%s", name);
    if (wi){
        fputc (' ', log_stream);
        fprintf (log_stream, "%d", logging_rank);
    }
    fprintf (log_stream, "] %*s", lindent, "");

    if (priority == SC_LP_TRACE) {
        char                bn[BUFSIZ], *bp;

        snprintf (bn, BUFSIZ, "%s", filename);
        bp = basename (bn);
        fprintf (log_stream, "%s:%d ", bp, lineno);
    }

    fputs (msg, log_stream);
    fflush (log_stream);
}

static void
sc_log_handler (FILE * log_stream, const char *filename, int lineno,
                int package, int category, int priority, const char *msg)
{
    log_handler("libsc",log_stream,filename,lineno,package,category,priority,msg);
}

static void
p4est_log_handler (FILE * log_stream, const char *filename, int lineno,
                int package, int category, int priority, const char *msg)
{
    log_handler("p4est",log_stream,filename,lineno,package,category,priority,msg);
}

static void
fclaw_log_handler (FILE * log_stream, const char *filename, int lineno,
                int package, int category, int priority, const char *msg)
{
    log_handler("fclaw",log_stream,filename,lineno,package,category,priority,msg);
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
    fclaw_global_productionf ("%-*s %s\n", w, "F77", FCLAW_F77);
    fclaw_global_productionf ("%-*s %s\n", w, "FFLAGS", FCLAW_FFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "CC", FCLAW_CC);
    fclaw_global_productionf ("%-*s %s\n", w, "CFLAGS", FCLAW_CFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "CXX", FCLAW_CXX);
    fclaw_global_productionf ("%-*s %s\n", w, "CXXFLAGS", FCLAW_CXXFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "LDFLAGS", FCLAW_LDFLAGS);
    fclaw_global_productionf ("%-*s %s\n", w, "FLIBS", FCLAW_FLIBS);
    fclaw_global_productionf ("%-*s %s\n", w, "LIBS", FCLAW_LIBS);
}

fclaw_app_t *
fclaw_app_new_on_comm (sc_MPI_Comm mpicomm, int *argc, char ***argv, void *user)
{
    //TODO seperate intialize from creating new app (makes testing difficult)
#ifdef FCLAW_ENABLE_DEBUG
    const int LP_lib = SC_LP_INFO;
    const int LP_fclaw = SC_LP_DEBUG;
#else
    const int LP_lib = SC_LP_ESSENTIAL;
    const int LP_fclaw = SC_LP_PRODUCTION;
#endif
    int mpiret;
    fclaw_app_t *a;

    sc_init (mpicomm, 1, 1, sc_log_handler, LP_lib);
    p4est_init (p4est_log_handler, LP_lib);
    fclaw_init (fclaw_log_handler, LP_fclaw);

    a = FCLAW_ALLOC (fclaw_app_t, 1);
    a->mpicomm = mpicomm;
    mpiret = sc_MPI_Comm_size (a->mpicomm, &a->mpisize);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Comm_rank (a->mpicomm, &a->mpirank);
    SC_CHECK_MPI (mpiret);
    logging_rank = a->mpirank;

    srand (a->mpirank);
    a->first_arg = -1;
    a->argc = argc;
    a->argv = argv;
    a->user = user;
    a->opt = sc_options_new ((*argv)[0]);
    sc_options_set_spacing (a->opt, 40, 56);
    a->core_registered = 0;

    a->opt_pkg = sc_array_new (sizeof (fclaw_app_options_t));

    a->configdir = fclaw_configdir;
    a->env_configdir = fclaw_env_configdir;

    a->attributes = sc_keyvalue_new ();

    return a;
}

fclaw_app_t *
fclaw_app_new (int *argc, char ***argv, void *user)
{
    int mpiret;
    sc_MPI_Comm mpicomm;

    mpicomm = sc_MPI_COMM_WORLD;

    mpiret = sc_MPI_Init (argc, argv);
    SC_CHECK_MPI (mpiret);

    return fclaw_app_new_on_comm(mpicomm, argc, argv, user);
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

    /* destroy central structures */
    sc_keyvalue_destroy (a->attributes);
    sc_options_destroy (a->opt);

    /* let the options packages clean up their memory */
    for (zz = a->opt_pkg->elem_count; zz > 0; --zz)
    {
        ao = (fclaw_app_options_t *) sc_array_index (a->opt_pkg, zz - 1);
        FCLAW_ASSERT (ao != NULL);
        if (ao->vt.options_destroy != NULL)
        {
            ao->vt.options_destroy (a, ao->package, ao->registered);
        }
        FCLAW_FREE (ao->section);
        FCLAW_FREE (ao->configfile);
    }
    sc_array_destroy (a->opt_pkg);

    FCLAW_FREE (a);

    sc_finalize ();

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);
}

void *
fclaw_app_set_attribute (fclaw_app_t * a, const char *name, void *attribute)
{
    void *previous;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (a->attributes != NULL);
    FCLAW_ASSERT (name != NULL);

    /* This is ugly and twice as expensive as a good solution.
     * Currently the sc_keyvalue API is not optimized in this regard. */
    previous = sc_keyvalue_get_pointer (a->attributes, name, NULL);
    sc_keyvalue_set_pointer (a->attributes, name, attribute);

    return previous;
}

void *
fclaw_app_get_attribute (fclaw_app_t * a,
                         const char *name, void *default_return)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (a->attributes != NULL);
    FCLAW_ASSERT (name != NULL);

    return sc_keyvalue_get_pointer (a->attributes, name, default_return);
}

void
fclaw_app_set_configdir (fclaw_app_t * a, const char *configdir)
{
    FCLAW_ASSERT (a != NULL);

    a->configdir = configdir;
}

void
fclaw_app_set_env_configdir (fclaw_app_t * a, const char *env_configdir)
{
    FCLAW_ASSERT (a != NULL);

    a->env_configdir = env_configdir;
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

/** This is the internal state of an options structure for core variables. */
typedef struct fclaw_options_core
{
    int print_help;        /**< Option variable to activate help message */
    int print_version;     /**< Option variable to print the version */
    int print_options;     /**< Opiton variable to print options and exit */
    int fclaw_verbosity;   /**< Option variable for ForestClaw verbosity */
    int lib_verbosity;     /**< Option variable for p4est, sc, and others */
    sc_keyvalue_t *kv_verbosity;      /**< Holds key-values for log levels */

    /* this is just for ForestClaw debugging, no need to adopt elsewhere */
    int is_registered;     /**< Internal variable to double-check the flow */
}
fclaw_options_core_t;

static void *
options_register_core (fclaw_app_t * a, void *package, sc_options_t * opt)
{
    sc_keyvalue_t *kv;
    fclaw_options_core_t *core = (fclaw_options_core_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    /* allocated storage for this package's option values */
    FCLAW_ASSERT (core != NULL);
    FCLAW_ASSERT (!core->is_registered);

    /* this key-value pair understands the verbosity levels */
    kv = core->kv_verbosity = sc_keyvalue_new ();
    sc_keyvalue_set_int (kv, "default", FCLAW_VERBOSITY_DEFAULT);
    sc_keyvalue_set_int (kv, "debug", FCLAW_VERBOSITY_DEBUG);
    sc_keyvalue_set_int (kv, "info", FCLAW_VERBOSITY_INFO);
    sc_keyvalue_set_int (kv, "production", FCLAW_VERBOSITY_PRODUCTION);
    sc_keyvalue_set_int (kv, "essential", FCLAW_VERBOSITY_ESSENTIAL);
    sc_keyvalue_set_int (kv, "error", FCLAW_VERBOSITY_ERROR);
    sc_keyvalue_set_int (kv, "silent", FCLAW_VERBOSITY_SILENT);

    /* set the options for the core package */
    sc_options_add_switch (opt, 'h', "help", &core->print_help,
                           "Print usage information");
    sc_options_add_switch (opt, 'v', "version", &core->print_version,
                           "Print ForestClaw version");
    sc_options_add_switch (opt, 'P', "print-options", &core->print_options,
                           "Print option values and exit");
    sc_options_add_keyvalue (opt, 'V', "verbosity", &core->fclaw_verbosity,
                             "default", kv, "Set ForestClaw verbosity");
    sc_options_add_keyvalue (opt, '\0', "lib-verbosity", &core->lib_verbosity,
                             "essential", kv, "Set verbosity for libraries");
    sc_options_add_inifile (opt, 'F', "configfile",
                            "Optional configuration file");

    /* we do not need to work with the return value */
    core->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
options_postprocess_core (fclaw_app_t * a, void *package, void *registered)
{
    fclaw_options_core_t *core = (fclaw_options_core_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */

    /* postprocess this package */
    FCLAW_ASSERT (core != NULL);
    FCLAW_ASSERT (core->is_registered);

    /* go through this packages options */
    sc_package_set_verbosity (sc_package_id, core->lib_verbosity);
    sc_package_set_verbosity (p4est_package_id, core->lib_verbosity);
    sc_package_set_verbosity (fclaw_get_package_id (), core->fclaw_verbosity);

    /* print help and/or version information and exit gracefully */
    if (core->print_version)
    {
        fclaw_global_essentialf ("ForestClaw version %s\n",
                                 FCLAW_PACKAGE_VERSION);
        return FCLAW_EXIT_QUIET;
    }
    if (core->print_options)
    {
        return FCLAW_EXIT_PRINT;
    }
    if (core->print_help)
    {
        return FCLAW_EXIT_USAGE;
    }

    /* at this point there are no errors to report */
    return FCLAW_NOEXIT;
}

static void
options_destroy_core (fclaw_app_t * a, void *package, void *registered)
{
    fclaw_options_core_t *core = (fclaw_options_core_t *) package;

    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* free this package */
    FCLAW_ASSERT (core != NULL);
    FCLAW_ASSERT (core->is_registered);
    FCLAW_ASSERT (core->kv_verbosity != NULL);
    sc_keyvalue_destroy (core->kv_verbosity);

    FCLAW_FREE (core);
}

static const fclaw_app_options_vtable_t options_vtable_core = {
    options_register_core,
    options_postprocess_core,
    NULL,
    options_destroy_core
};

int
fclaw_app_options_core_registered (fclaw_app_t * a)
{
    return a->core_registered;
}

void
fclaw_app_options_register_core (fclaw_app_t * a, const char *configfile)
{
    fclaw_options_core_t *core;

    FCLAW_ASSERT (a != NULL);

    /* allocate storage for core's option values */
    /* we will free it in the options_destroy callback */
    core = FCLAW_ALLOC_ZERO (fclaw_options_core_t, 1);

    /* sneaking the version string into the package pointer */
    /* when there are more parameters to pass, create a structure to pass */
    fclaw_app_options_register (a, NULL, configfile, &options_vtable_core,
                                core);
    a->core_registered = 1;
}


void fclaw_app_print_options(fclaw_app_t *app)
{
        sc_options_print_summary (fclaw_get_package_id (),
                                  FCLAW_VERBOSITY_ESSENTIAL, app->opt);    
}

fclaw_exit_type_t
fclaw_app_options_parse (fclaw_app_t * a, int *first_arg,
                         const char *savefile)
{
    int retval;
    size_t zz;
    fclaw_exit_type_t vexit;
    fclaw_app_options_t *ao;

    FCLAW_ASSERT (a != NULL);

    /* TODO: read configuration files */

    /* parse command line options with given priority for errors */
    a->first_arg =
        sc_options_parse (fclaw_get_package_id (), FCLAW_VERBOSITY_ERROR,
                          a->opt, *a->argc, *a->argv);

    /* check for option and parameter errors */
    if (a->first_arg < 0)
    {
        /* option processing was not successful */
        vexit = FCLAW_EXIT_ERROR;
    }
    else
    {
        /* go through options packages for further processing and verification */
        vexit = FCLAW_NOEXIT;
        for (zz = 0; zz < a->opt_pkg->elem_count; ++zz)
        {
            fclaw_exit_type_t aoexit;

            ao = (fclaw_app_options_t *) sc_array_index (a->opt_pkg, zz);
            FCLAW_ASSERT (ao != NULL);
            if (ao->vt.options_postprocess != NULL)
            {
                aoexit = ao->vt.options_postprocess (a, ao->package,
                                                     ao->registered);
                vexit = SC_MAX (aoexit, vexit);
            }
            if (ao->vt.options_check != NULL)
            {
                aoexit = ao->vt.options_check (a, ao->package,
                                               ao->registered);
                vexit = SC_MAX (aoexit, vexit);
            }
        }
    }

    /* let's see what we print */
    /* rationale: only use ESSENTIAL for the primary purpose of an exit condition
     *            only use PRODUCTION for really useful information
     *            partially redundant output can go with INFO
     */
    switch (vexit)
    {
    case FCLAW_NOEXIT:
        fclaw_global_infof ("Option parsing successful\n");
        /* Ok maybe this is essential for now.
         * Note: standard log level for a run used for producing results
         *       should be PRODUCTION.  ESSENTIAL is only for those who
         *       choose to ignore potentially important information.
         */
        sc_options_print_summary (fclaw_get_package_id (),
                                  FCLAW_VERBOSITY_ESSENTIAL, a->opt);
        break;
    case FCLAW_EXIT_QUIET:
        /* we assume that the application has printed or will print something */
        break;
    case FCLAW_EXIT_PRINT:
        /* it has been requested to print values of all options and then exit */
        sc_options_print_summary (fclaw_get_package_id (),
                                  FCLAW_VERBOSITY_ESSENTIAL, a->opt);
        fclaw_global_infof ("Terminating after printing option values\n");
        break;
    case FCLAW_EXIT_USAGE:
        /* we assume that the application has printed or will print something */
        /* but it has been specifically requested to print usage information */
        sc_options_print_usage (fclaw_get_package_id (),
                                FCLAW_VERBOSITY_ESSENTIAL, a->opt, NULL);
        fclaw_global_infof ("Terminating program by request\n");
        break;
    case FCLAW_EXIT_ERROR:
        /* some error has been encountered */
        fclaw_global_infof ("Configuration / option parsing failed\n");
        sc_options_print_usage (fclaw_get_package_id (),
                                FCLAW_VERBOSITY_PRODUCTION, a->opt, NULL);
        fclaw_global_errorf ("Terminating program on option error\n");
        break;
    default:
        SC_ABORT_NOT_REACHED ();
    }

    /* print configuration if so desired */
    if (vexit != FCLAW_EXIT_ERROR && sc_is_root () && savefile != NULL)
    {
        retval = sc_options_save (fclaw_get_package_id (),
                                  FCLAW_VERBOSITY_ERROR, a->opt, savefile);
        if (retval)
        {
            vexit = FCLAW_EXIT_ERROR;
            fclaw_global_infof ("Unable to save options to \"%s\"\n",
                                savefile);
        }
    }

    /* we are done */
    if (first_arg != NULL)
    {
        *first_arg = a->first_arg;
    }
    return vexit;
}

sc_MPI_Comm
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

static fclaw_pointer_map_t* packing_vtables = NULL;

void fclaw_app_register_options_packing_vtable(const char*name,fclaw_packing_vtable_t* vtable){
    if(packing_vtables == NULL)
    {
        packing_vtables = fclaw_pointer_map_new();
    }
    fclaw_pointer_map_insert(packing_vtables, name, vtable, NULL);
}

fclaw_packing_vtable_t* fclaw_app_get_options_packing_vtable(const char*name){
    if(packing_vtables == NULL)
    {
        packing_vtables = fclaw_pointer_map_new();
    }
    return (fclaw_packing_vtable_t*) fclaw_pointer_map_get(packing_vtables,name);
}

/*** which of the following do we need? ***/

#if 0

void
fclaw2d_global_log (int log_priority, const char *message)
{
    /* TODO: establish an fclaw_package_id */
    SC_GEN_LOG (sc_package_id, SC_LC_GLOBAL, log_priority, message);
}

void *
fclaw2d_alloc (size_t size)
{
    return sc_malloc (p4est_package_id, size);
}

void *
fclaw2d_calloc (size_t nmemb, size_t size)
{
    return sc_calloc (p4est_package_id, nmemb, size);
}

void *
fclaw2d_realloc (void *ptr, size_t size)
{
    return sc_realloc (p4est_package_id, ptr, size);
}

void
fclaw2d_free (void *ptr)
{
    sc_free (p4est_package_id, ptr);
}

#endif /* 0 */
