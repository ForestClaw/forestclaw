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
 *    or unexpected error conditions.  Conversely, log messages that
 *    help in debugging and clutter the output a lot can use
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
 *     FCLAW_VERBOSITY_DEFAULT, DEBUG, INFO, PRODUCTION, ESSENTIAL, SILENT
 * (although I suppose we could use others; see below for the list of SC_LP_ values).
 *
 * * PRODUCTION messages should log a few lines for a major api function.
 *              A good example would be the numerical error at the
 *              end of the program, timings, printing important parameters.
 * * ESSENTIAL messages must only cause a few lines for the whole run, if at all.
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
    FCLAW_VERBOSITY_ESSENTIAL = SC_LP_ESSENTIAL,
    FCLAW_VERBOSITY_PRODUCTION = SC_LP_PRODUCTION,
    FCLAW_VERBOSITY_INFO = SC_LP_INFO,
    FCLAW_VERBOSITY_DEBUG = SC_LP_DEBUG,
    FCLAW_VERBOSITY_DEFAULT
}
fclaw_verbosity_t;

/** An application container whose use is optional. */
typedef struct fclaw_app fclaw_app_t;

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
 * \param [in] priority     The log level can be FCLAW_VERBOSITY_ESSENTIAL,
 *                          _PRODUCTION, _INFO, and _DEBUG, (see SC_LP_*).
 * \param [in] fmt          Format string as in printf.
 */
void fclaw_logf (int category, int priority, const char *fmt, ...)
#ifndef FCLAW_DOXYGEN
    __attribute__ ((format (printf, 3, 4)))
#endif
    ;
/** Print a message only on the root processor with priority FCLAW_VERBOSITY_ESSENTIAL.
 * This priority must be used \a very sparingly, such as for a version number,
 * package registration, and error messages.
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
 * This can be influenced at compile time with --enable-logging=SC_LP_DEBUG
 * for example, but this is somewhat clumsy and usually unnecessary since this
 * option does not differentiate between the forestclaw and its submodules.
 * It is possible and encouraged to change the levels with sc_package_set_verbosity.
 * Attempts to reduce them (i.e., to cause more verbosity) at runtime are ignored.
 * \param [in,out] argc         Command line argument count.
 * \param [in,out] argv         Command line arguments.
 * \param [in,out] user         Will not be changed by our code.
 * \return            An allocated and initialized application object.
 */
fclaw_app_t *fclaw_app_new (int *argc, char ***argv, void *user);

/** Close down all systems that were setup in fclaw_init.
 * If a keyvalue structure has been added to a->opt, it is destroyed too.
 */
void fclaw_app_destroy (fclaw_app_t * a);

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
MPI_Comm fclaw_app_get_mpi_size_rank (fclaw_app_t * a,
                                      int *mpisize, int *mpirank);

/** Return a pointer to the options structure.
 * TODO: aiming to provide an encapsulation that will not need this function.
 */
sc_options_t *fclaw_app_get_options (fclaw_app_t * a);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW_BASE_H */
