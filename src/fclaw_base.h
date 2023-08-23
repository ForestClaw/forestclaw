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

/** \file
 *
 * Basic include directives and declarations for ForestClaw.
 *
 * This file helps to integrate with p4est and libsc and provides general
 * macros.  It should be included regardless of the space dimension that
 * the code is compiled for.
 *
 * This file also provides logging and memory allocation functions, and a
 * convenience application object that should reduce boilerplate code.
 */

#ifndef FCLAW_BASE_H
#define FCLAW_BASE_H

/* include config headers */

#include <fclaw_config.h>
#include <sc_config.h>
#if \
  (defined (FCLAW_ENABLE_MPI) && !defined (SC_ENABLE_MPI)) || \
  (!defined (FCLAW_ENABLE_MPI) && defined (SC_ENABLE_MPI))
#error "MPI configured differently in ForestClaw and libsc"
#endif
#if \
  (defined (FCLAW_ENABLE_MPIIO) && !defined (SC_ENABLE_MPIIO)) || \
  (!defined (FCLAW_ENABLE_MPIIO) && defined (SC_ENABLE_MPIIO))
#error "MPI I/O configured differently in ForestClaw and libsc"
#endif
#include <p4est_base.h>
#define _fclaw_const _sc_const          /**< Replaces the const keyword. */
#define _fclaw_restrict _sc_restrict    /**< Replaces the restrict keyword. */

/* include some standard headers */

#include <sc_options.h>

/* define F77 and FC name-mangling if autoconf fails to do so */
#ifndef FCLAW_F77_FUNC
#define FCLAW_F77_FUNC(name,NAME) name ## _
#endif

#ifndef FCLAW_F77_FUNC_
#define FCLAW_F77_FUNC_(name,NAME) name ## _
#endif

#ifndef FCLAW_FC_FUNC
#define FCLAW_FC_FUNC(name,NAME) name ## _
#endif

#ifndef FCLAW_FC_FUNC_
#define FCLAW_FC_FUNC_(name,NAME) name ## _
#endif

/* start declarations */

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/** These are our log levels, sorted by priority.
 * These levels have two different roles to play:
 *
 * 1. Each logging function is suffixed with its priority.
 *    This means that FCLAW_VERBOSITY_ESSENTIAL should only be used
 *    for very important messages, such as a version number on startup,
 *    and FCLAW_VERBOSITY_ERROR for unexpected error conditions.  Conversely,
 *    log messages that help in debugging and clutter the output a lot can use
 *    the FCLAW_VERBOSITY_DEBUG priority.
 * 2. The function \b sc_package_set_verbosity
 *    accepts a priority level below which all messages will be ignored.
 *    Eventually, the \ref fclaw_app_t object will read the desired
 *    verbosity from the command line and/or configuration files.
 *    As an example, if the verbosity is set to FCLAW_VERBOSITY_INFO, then
 *    debug messages will not be printed, but essential messages will.
 *
 * The verbosity of libraries, like p4est and sc, can be controlled
 * individually, and separately from the verbosity of forestclaw.
 * For production runs, try to set the verbosity to PRODUCTION or ESSENTIAL.
 *
 * Valid ForestClaw settings:
 *     FCLAW_VERBOSITY_DEFAULT, DEBUG, INFO, PRODUCTION, ESSENTIAL, ERROR, SILENT
 * (although I suppose we could use others; see below for the list of SC_LP_ values).
 *
 * * PRODUCTION messages should log a few lines for a major api function.
 *              A good example would be the numerical error at the
 *              end of the program, timings, printing important parameters.
 * * ESSENTIAL messages must only cause a few lines for the whole run, if at all.
 * * ERROR messages should only be printed when something goes wrong.  Setting the
 *              verbosity to ERROR omits even essential messages and prints errors only.
 * * DEFAULT depends on the configure options `--enable-debug` and `--enable_logging=...`
 *              With `--enable-debug`, this will be a very chatty at level DEBUG.
 *              Otherwise it is INFO, if not overridden with `--enable-logging=...`
 *
 * This is copied from \b sc.h for reference:
 *
 * For debugging, try INFO or PRODUCTION for libraries and DEBUG for forestclaw.
 *
 *     #define SC_LP_DEFAULT   (-1)    this selects the SC default threshold
 *     #define SC_LP_ALWAYS      0     this will log everything
 *     #define SC_LP_TRACE       1     this will prefix file and line number
 *     #define SC_LP_DEBUG       2     any information on the internal state
 *     #define SC_LP_VERBOSE     3     information on conditions, decisions
 *     #define SC_LP_INFO        4     the main things a function is doing
 *     #define SC_LP_STATISTICS  5     important for consistency/performance
 *     #define SC_LP_PRODUCTION  6     a few lines for a major api function
 *     #define SC_LP_ESSENTIAL   7     this logs a few lines max per program
 *     #define SC_LP_ERROR       8     this logs errors only
 *     #define SC_LP_SILENT      9     this never logs anything
 */
typedef enum fclaw_verbosity
{
    FCLAW_VERBOSITY_SILENT = SC_LP_SILENT,
    FCLAW_VERBOSITY_ERROR = SC_LP_ERROR,
    FCLAW_VERBOSITY_ESSENTIAL = SC_LP_ESSENTIAL,
    FCLAW_VERBOSITY_PRODUCTION = SC_LP_PRODUCTION,
    FCLAW_VERBOSITY_INFO = SC_LP_INFO,
    FCLAW_VERBOSITY_DEBUG = SC_LP_DEBUG,
    FCLAW_VERBOSITY_DEFAULT
}
fclaw_verbosity_t;

/** This enumeration defines values to be returned from option checking.
 * The order of these constants is important; we rely on it internally.
 * */
typedef enum fclaw_exit_type
{
    FCLAW_NOEXIT,               /**< We are completely clean and may continue.
                                     By the C standard, this constant has value 0. */
    FCLAW_EXIT_QUIET,           /**< We decided to to terminate quietly. */
    FCLAW_EXIT_PRINT,           /**< We will print all option values and then quit. */
    FCLAW_EXIT_USAGE,           /**< We are to do an informative job and then quit.
                                     ForestClaw will print a usage message. */
    FCLAW_EXIT_ERROR            /**< We have encountered an error and quit noisily. */
}
fclaw_exit_type_t;

/** Making it more explicit when we are using boolean values. */
typedef int int_t;

/** This function turns an exit type into a value suitable for ending main ().
 * \param [in] vexit    Exit type from enumeration.
 * \return              A value suitable for returning from main () or
 *                      as an argument to exit ().
 */
int fclaw_app_exit_type_to_status (fclaw_exit_type_t vexit);

/** An application container whose use is optional. */
typedef struct fclaw_app fclaw_app_t;

/** Callback function type for registering options per-package. */
typedef void *(*fclaw_app_options_register_t) (fclaw_app_t * a,
                                               void *package,
                                               sc_options_t * options);

/** Callback function type for postprocessing options per-package.
 * This function as called after option parsing.
 * Its return value can be used to report an error in option processing. */
typedef fclaw_exit_type_t
    (*fclaw_app_options_postprocess_t) (fclaw_app_t * a,
                                        void *package, void *registered);

/** Callback function type for checking for option errors per-package.
 * This function as called after option parsing and postprocessing.
 * Its return value is used to report an error in option processing. */
typedef fclaw_exit_type_t
    (*fclaw_app_options_check_t) (fclaw_app_t * a,
                                  void *package, void *registered);

/** Virtual function to destroy any memory created on registration.
 * This is called at the end of the program and reverses any allocation
 * that has been done on option registration. */
typedef void (*fclaw_app_options_destroy_t) (fclaw_app_t * a, void *package,
                                             void *registered);

/** An extended table of virtual functions for package-wise option parsing. */
typedef struct fclaw_app_options_vtable
{
    fclaw_app_options_register_t options_register;       /**< Calls sc_options_add_*. */
    fclaw_app_options_postprocess_t options_postprocess; /**< Called after parsing. */
    fclaw_app_options_check_t options_check;     /**< Reports on validity of option values. */
    fclaw_app_options_destroy_t options_destroy;     /**< Called at end of program. */
}
fclaw_app_options_vtable_t;

/* macros for memory allocation, will abort if out of memory */

/** Allocate memory for n objects of a certain type. */
#define FCLAW_ALLOC(t,n)          (t *) sc_malloc (fclaw_get_package_id (), \
                                                   (n) * sizeof(t))
/** Allocate memory for n objects of a certain type and initialize to all zeros. */
#define FCLAW_ALLOC_ZERO(t,n)     (t *) sc_calloc (fclaw_get_package_id (), \
                                                   (size_t) (n), sizeof(t))
/** Reallocate memory of a certain type to n objects. */
#define FCLAW_REALLOC(p,t,n)      (t *) sc_realloc (fclaw_get_package_id (), \
                                                    (p), (n) * sizeof(t))
/** Allocate memory for a given string and copy the string into it. */
#define FCLAW_STRDUP(s)                 sc_strdup (fclaw_get_package_id (), (s))
/** Free memory from any of the FCLAW_ allocation functions.
 * This macro must not be used to free memory from p4est, sc, or other packages. */
#define FCLAW_FREE(p)                   sc_free (fclaw_get_package_id (), (p))

/* some error checking */
#ifdef FCLAW_ENABLE_DEBUG
#define FCLAW_ASSERT(c) SC_CHECK_ABORT ((c), "Assertion '" #c "'")
#define FCLAW_EXECUTE_ASSERT_FALSE(expression)                          \
  do { int _fclaw_i = (int) (expression);                               \
       SC_CHECK_ABORT (!_fclaw_i, "Expected false: '" #expression "'"); \
  } while (0)
#define FCLAW_EXECUTE_ASSERT_TRUE(expression)                           \
  do { int _fclaw_i = (int) (expression);                               \
       SC_CHECK_ABORT (_fclaw_i, "Expected true: '" #expression "'");   \
  } while (0)
#else
/** Without --enable-debug, assertions do nothing. */
#define FCLAW_ASSERT(c) SC_NOOP ()
/** Calculate a given expression.  No side effects. */
#define FCLAW_EXECUTE_ASSERT_FALSE(expression) \
  do { (void) (expression); } while (0)
/** Calculate a given expression.  No side effects. */
#define FCLAW_EXECUTE_ASSERT_TRUE(expression) \
  do { (void) (expression); } while (0)
#endif
#define FCLAW_ABORT(error_message)  SC_ABORT (error_message)

/** Query the package identity for use with libsc functions.
 * This is -1 until fclaw_init or fclaw_app_new has been called.
 * \return              A package identifier to pass to log functions etc.
 */
int fclaw_get_package_id (void);

/* *INDENT-OFF* */
/** Function for printing log messages that rely on sc_logv (see \b sc.h).
 * \param [in] category     The category SC_LC_GLOBAL means that only rank 0
 *                          prints any output.  SC_LC_NORMAL means that all
 *                          ranks prints output, which is generally undesired
 *                          in production runs (that is, without
 *                          --enable-debug).
 * \param [in] priority     The log level can be FCLAW_VERBOSITY_ERROR, _ESSENTIAL,
 *                          _PRODUCTION, _INFO, and _DEBUG, (see SC_LP_*).
 * \param [in] fmt          Format string as in printf.
 */
void fclaw_logf (int category, int priority, const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 3, 4)))
#endif
    ;
/** Print a message only on the root processor with priority FCLAW_VERBOSITY_ERROR.
 * This priority must only used for unexpected error conditions that threaten
 * the successful completion of the program.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_global_errorf (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/** Print a message on the calling processor with priority FCLAW_VERBOSITY_ERROR.
 * For errors that likely occur similarly on all processors, please use \ref
 * fclaw_global_errorf instead to avoid flooding the output.  This function
 * must only used for unexpected individual-processor error conditions that
 * threaten the successful completion of the program.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_errorf (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/** Print a message only on the root processor with priority FCLAW_VERBOSITY_ESSENTIAL.
 * This priority must be used \a very sparingly, such as for a version number
 * or a confirmation of package registration.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_global_essentialf (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/** Print a message only on the root processor with priority FCLAW_VERBOSITY_PRODUCTION.
 * This priority is useful for high-level results or statistics output.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_global_productionf (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/** Print a message only on the root processor with priority FCLAW_VERBOSITY_INFO.
 * This priority is useful for less important message on the flow of the program.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_global_infof (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/** Print a processor-local message with priority FCLAW_VERBOSITY_INFO.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_infof (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/** Print a processor-local message with priority FCLAW_VERBOSITY_DEBUG.
 * \param [in] fmt      A printf-style format string.
 */
void fclaw_debugf (const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 1, 2)))
#endif
    ;
/* *INDENT-ON* */

/** Register ForestClaw with libsc and print version and variable information.
 * It is not necessary to call this function, but it makes the log output more
 * useful and separates ForestClaw's memory allocation from other packages.
 * This function is usually not called directly since it is automatically
 * called from \ref fclaw_app_new.
 * \param [in] log_handler   Declared in sc.h.  Usually, NULL is fine.
 * \param [in] log_threshold Declared in sc.h.  SC_LP_DEFAULT is fine.
 *                           You can also choose from log levels SC_LP_*.
 */
void fclaw_init (sc_log_handler_t log_handler, int log_threshold);

/** Call the (optional) init functions for sc, p4est, and ForestClaw.  The log
 * level is set as well and depends on the configure option `--enable-debug`.
 * With `--enable-debug`: DEBUG for ForestClaw, INFO for sc and p4est.  Without
 * `--enable-debug`: PRODUCTION for ForestClaw, ESSENTIAL for sc and p4est.
 * This can be influenced at compile time with `--enable-logging=SC_LP_DEBUG`
 * for example, but this is somewhat clumsy and usually unnecessary since this
 * option does not differentiate between the forestclaw and its submodules.
 * It is possible and encouraged to change the levels with \b sc_package_set_verbosity.
 * Attempts to reduce them (i.e., to cause more verbosity) at runtime are ignored.
 * \param [in,out] argc         Command line argument count.
 * \param [in,out] argv         Command line arguments.
 * \param [in,out] user         Will not be changed by our code.
 * \return            An allocated and initialized application object.
 */
fclaw_app_t *fclaw_app_new (int *argc, char ***argv, void *user);

/** Use when MPI is already intialized.
 * Call the (optional) init functions for sc, p4est, and ForestClaw.  The log
 * level is set as well and depends on the configure option `--enable-debug`.
 * With `--enable-debug`: DEBUG for ForestClaw, INFO for sc and p4est.  Without
 * `--enable-debug`: PRODUCTION for ForestClaw, ESSENTIAL for sc and p4est.
 * This can be influenced at compile time with `--enable-logging=SC_LP_DEBUG`
 * for example, but this is somewhat clumsy and usually unnecessary since this
 * option does not differentiate between the forestclaw and its submodules.
 * It is possible and encouraged to change the levels with \b sc_package_set_verbosity.
 * Attempts to reduce them (i.e., to cause more verbosity) at runtime are ignored.
 * \param [in]     mpicomm      The MPI comm to initialize on.
 * \param [in,out] argc         Command line argument count.
 * \param [in,out] argv         Command line arguments.
 * \param [in,out] user         Will not be changed by our code.
 * \return            An allocated and initialized application object.
 */
fclaw_app_t *fclaw_app_new_on_comm (sc_MPI_Comm mpicomm, int *argc, char ***argv, void *user);


/** Close down all systems that were setup in fclaw_init.
 * If a keyvalue structure has been added to a->opt, it is destroyed too.
 * \param [in,out] a            This application is cleaned up.
 */
void fclaw_app_destroy (fclaw_app_t * a);

/** Set or replace a named attribute in an application.
 * If there is already an attribute stored under the given name,
 * return the old content and overwrite it.
 * The application object does no interpretation of the content.
 * The content can be accessed again with \ref fclaw_app_get_attribute.
 * All attributes are removed automatically by \ref fclaw_app_destroy.
 * \param [in,out] a            Valid application object.
 * \param [in] name             The name, or lookup key, for an attribute.
 *                              It must not go out of scope while in use.
 * \param [in] attribute        The new content, or value, of an attribute.
 *                              Setting this to NULL does save NULL as value
 *                              and the attribute still continues to exist.
 * \return                      The value previously stored under this name,
 *                              or NULL if the attribute does not exist yet.
 */
void *fclaw_app_set_attribute (fclaw_app_t * a,
                               const char *name, void *attribute);

/** Retrieve a named attribute stored in an application.
 * If the attribute does not exist, return the specified default value.
 * \param [in] a                Valid application object.
 * \param [in] name             The name, or lookup key, for an attribute.
 *                              It must not go out of scope while in use.
 * \param [in] default_return   This value is return if name does not exist.
 * \return                      The value previously stored under this name,
 *                              or \b default_return if it does not exist.
 */
void *fclaw_app_get_attribute (fclaw_app_t * a,
                               const char *name, void *default_return);

/** Change the user's configuration directory.
 * This is understood relative to the user's home directory.
 * After \ref fclaw_app_new, it is set to ".forestclaw".
 * \param [in,out] a            Valid application object.
 * \param [in] configdir        New value for the per-user configuration directory.
 *                              This string is passed, not copied, so it must not
 *                              go out of scope.  String constants are fine.
 *                              It may be NULL to deactivate the feature.
 */
void fclaw_app_set_configdir (fclaw_app_t * a, const char *configdir);

/** Change the environment variable used to find another configuration directory.
 * After \ref fclaw_app_new, it is set to "FCLAW_INI_DIR".
 * \param [in,out] a            Valid application object.
 * \param [in] env_configdir    New value for the environment variable.
 *                              This string is passed, not copied, so it must not
 *                              go out of scope.  String constants are fine.
 *                              It may be NULL to deactivate the feature.
 */
void fclaw_app_set_env_configdir (fclaw_app_t * a, const char *env_configdir);

/** Register an options package with an application.
 * This function calls the virtual function for registering its options.
 * \param [in,out] a            A valid application object from \ref fclaw_app_new.
 *                              The callback functions from \b vt may access its
 *                              user data, or work with \b package, or both.
 * \param [in] section          This variable decides about the options being
 *                              under the default "Options" section in a possible
 *                              configuration file, or under its own section.
 *                              If this is NULL, the options will be used with
 *                              the names used to add them to the options structure.
 *                              If not NULL, we will use \b sc_options_add_suboptions:
 *                              The true option name will be prefixed with \b section,
 *                              and appear under [section] in a .ini-style file.
 *                              Furthermore, its short options will be disabled
 *                              and long options prefixed with section:.
 * \param [in] configfile       IF not NULL, the name of a configuration file without
 *                              the path (but with the ending).  The file is read before
 *                              option parsing occurs, so the command line overrides.
 *                              TODO: this feature is not yet active.
 * \param [in] vt               Functions for options processing.  At least the
 *                              member \b options_register must be non-NULL.
 *                              If any of the virtual function produce console output,
 *                              they shall use the fclaw_global_* logging functions,
 *                              or custom functions that only print on rank zero.
 * \param [in] package          Any kind of context, or NULL, as a convenience for
 *                              the caller.  Is passed to all options callbacks in
 *                              \b vt.  This can be used in addition, or instead of,
 *                              accessing the user pointer from the application object.
 *                              The application may define that package is owned by
 *                              the internal options processing and free it in the
 *                              options_destroy callback, or impose other conventions.
 */
void fclaw_app_options_register (fclaw_app_t * a,
                                 const char *section, const char *configfile,
                                 const fclaw_app_options_vtable_t * vt,
                                 void *package);

/**
 * \brief Check if core options have been registered for this app
 * 
 * \param [in,out] a            A valid application object.
 * \return int                  1 if registered, 0 if not
 */
int
fclaw_app_options_core_registered (fclaw_app_t * a);

/** Register a central convenience options package with default behavior.
 * It is just an example and completely fine not to use this function.
 * This is not a replacement for calling \ref fclaw_app_options_register,
 * which may be used any number of times for other custom options packages.
 * It merely calls \ref fclaw_app_options_register with predefined operations.
 * This options package provides the following options:
 *
 *     -?, --help               Print a usages message for all options and exit.
 *     -v, --version            Print a version string and exit.
 *     -V, --verbosity=...      Set the verbosity for ForestClaw; a string in
 *                              \a lowercase letters without the prefix FCLAW_VERBOSITY_.
 *     --lib-verbosity=...      Like verbosity, but for the libraries p4est and sc.
 *     -F, --configfile=...     The name/path to a configuration file that is read
 *                              while parsing the options from the command line.
 *                              Handling is different from the configfile argument
 *                              to fclaw_app_options_register in that it must exist
 *                              if asked for, that no default paths are tried for
 *                              its location, and that it is read during option
 *                              parsing, possibly modifying values set only recently.
 * \param [in,out] a            A valid application object.
 * \param [in] configfile       If not NULL, an .ini-style configuration file is read
 *                              before option parsing.  This is its name without path.
 */
void fclaw_app_options_register_core (fclaw_app_t * a,
                                      const char *configfile);

/** Parse the command line options.
 * This function will loop through all registered packages and call the functions
 * for postprocessing and checking options, and will abort the program if any
 * error occurs.
 * TODO: This function shall read default configuration files before parsing.
 * \param [in] a         Initialized forestclaw application.
 * \param [out] first_arg       If not NULL, position of first non-option argument.
 * \param [in] savefile         If not NULL, write options to this file after parsing,
 *                              but only if the exit type is not an error.
 * \return               Whether to continue, exit gracefully, or exit with error.
 */
fclaw_exit_type_t fclaw_app_options_parse (fclaw_app_t * a, int *first_arg,
                                           const char *savefile);


/** Print out use options
 * \param [in] a         Initialized forestclaw application.
 * \return               This pointer is returned unchanged.
 */
void fclaw_app_print_options(fclaw_app_t *app);


/** Return the user pointer passed on \ref fclaw_app_new.
 * \param [in] a         Initialized forestclaw application.
 * \return               This pointer is returned unchanged.
 */
void *fclaw_app_get_user (fclaw_app_t * a);

/** Query MPI size and rank for an application.
 * \param [in] a        Initialized forestclaw application.
 * \param [out] mpisize Will be set to the size of the world communicator.
 * \param [out] mpirank Will be set to the size of the world communicator.
 * \return              The communicator that was used to setup \b a.
 */
sc_MPI_Comm fclaw_app_get_mpi_size_rank (fclaw_app_t * a,
                                         int *mpisize, int *mpirank);

/** Return a pointer to the options structure.
 * \deprecated TODO: We shall provide an interface that will not need this function.
 */
sc_options_t *fclaw_app_get_options (fclaw_app_t * a);

/**
 * @brief Set a logging prefix.
 *
 * This will prepend a [prefix] on all logging messages.
 * This is useful when running with two solvers.
 * 
 * @param prefix the logging prefix
 */
void fclaw_set_logging_prefix(const char* prefix);
#if 0

/*** rename the following names without the 2D ***/

///@}
/* ---------------------------------------------------------------------- */
///                      @name Fundamentals
/* ---------------------------------------------------------------------- */
///@{

/** Log a message only on rank zero. */
void fclaw2d_global_log (int log_priority, const char *message);

///@}
/* ---------------------------------------------------------------------- */
///                         @name Allocation
/* ---------------------------------------------------------------------- */
///@{

/**
 * @brief Allocate memory
 *
 * @param size the size in bytes
 * @return void* the newly allocated memory
 */
void *fclaw2d_alloc (size_t size);
/**
 * @brief Allocate a number of objects, intializes to 0
 *
 * @param nmemb the number of objects
 * @param size the size of an object (in bytes)
 * @return void* the newly allocated memory
 */
void *fclaw2d_calloc (size_t nmemb, size_t size);
/**
 * @brief Reallocate memory
 *
 * @param ptr the memory to be reallocated
 * @param size the size in bytes
 * @return void* the newly allocated memory
 */
void *fclaw2d_realloc (void *ptr, size_t size);
/**
 * @brief Free allocated memory
 *
 * @param ptr the memory to free
 */
void fclaw2d_free (void *ptr);
/**
 * @brief Allocate an array with given length
 *
 * @param t the type
 * @param n the length of the array
 * @return t* the newly allocateed memory
 */
#define FCLAW2D_ALLOC(t,n)      (t *) fclaw2d_alloc ((n) * sizeof (t))
/**
 * @brief Allocate an array with given length, initialized to zero
 *
 * @param t the type
 * @param n the length of the array
 * @return t* the newly allocateed memory
 */
#define FCLAW2D_ALLOC_ZERO(t,n) (t *) fclaw2d_calloc ((n), sizeof (t))
/**
 * @brief Allocate an array with given length, initialized to zero
 *
 * @param p pointer to memory to be reallocated
 * @param t the type
 * @param n the length of the array
 * @return t* the newly allocateed memory
 */
#define FCLAW2D_REALLOC(p,t,n)  (t *) fclaw2d_realloc ((p), (n) * sizeof (t))
/**
 * @brief Free allocated memory
 *
 * @param p the memory to free
 */
#define FCLAW2D_FREE(p)         fclaw2d_free (p)

#endif /* 0 */

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW_BASE_H */
