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
#include <forestclaw2d.h>
#include <fclaw_options.h>

int fclaw_options_read_from_file(sc_options_t* opt, int log_priority)
{
    int retval;
    retval = sc_options_load (sc_package_id, SC_LP_ALWAYS, opt,
                              "fclaw_options.ini");
    if (retval < 0)
    {
        fclaw2d_global_log (log_priority,      \
                            "Cannot find (or cannot open) fclaw_options.ini.\n");
    }
    else
    {
        fclaw2d_global_log (log_priority,      \
                            "Reading file fclaw_options.ini.\n");
    }
    return retval;
}


int fclaw_options_parse_command_line(sc_options_t * opt,
                                     int argc, char **argv,
                                     int log_priority)
{
    int retval;

    retval = sc_options_parse (sc_package_id, SC_LP_ERROR, opt, argc, argv);
    if (retval < 0)
    {
#if 0
        sc_options_print_usage (sc_package_id, log_priority, opt, NULL);
#endif
        fclaw2d_global_log (log_priority, "Command line option parsing failed.  " \
                            "Use --help option to see valid option settings.\n");
        return retval;
    }
    else
    {
        fclaw2d_global_log (log_priority, "Reading command line.\n");
    }
    if (sc_is_root ())
    {
        retval = sc_options_save (sc_package_id, SC_LP_ERROR, opt,
                                  "fclaw_options.ini.used");
        if (retval < 0)
        {
            fclaw2d_global_log (log_priority, "Cannot save fclaw_options.ini.used.  ");
            return retval;
        }
    }
    return retval;
}

void fclaw_options_print_summary(sc_options_t *opt, int log_priority)
{
    sc_options_print_summary (sc_package_id, log_priority, opt);
}

void
fclaw_options_add_int_array (sc_options_t * opt,
                             int opt_char, const char *opt_name,
                             const char **array_string,
                             const char *default_string,
                             int **int_array, int initial_length,
                             const char *help_string)
{
    *int_array = NULL;
    sc_options_add_string (opt, opt_char, opt_name,
                           array_string, default_string, help_string);
}

void
fclaw_options_add_double_array (sc_options_t * opt,
                                int opt_char, const char *opt_name,
                                const char **array_string,
                                const char *default_string,
                                double **double_array, int initial_length,
                                const char *help_string)
{
    *double_array = NULL;
    sc_options_add_string (opt, opt_char, opt_name,
                           array_string, default_string, help_string);
}

void
    fclaw_options_convert_int_array (const char *array_string,
                                     int **int_array, int new_length)
{
    int i;
    const char *beginptr;
    char *endptr;

    new_length = SC_MAX (new_length, 0);
    *int_array = FCLAW_REALLOC (*int_array, int, new_length);

    beginptr = array_string;
    for (i = 0; i < new_length; ++i)
    {
        if (beginptr == NULL)
        {
            (*int_array)[i] = 0;
        }
        else
        {
            (*int_array)[i] = (int) strtol (beginptr, &endptr, 10);
            beginptr = endptr;
        }
    }
}

void
    fclaw_options_convert_double_array (const char *array_string,
                                        double **double_array, int new_length)
{
    int i;
    const char *beginptr;
    char *endptr;

    new_length = SC_MAX (new_length, 0);
    *double_array = FCLAW_REALLOC (*double_array, double, new_length);

    beginptr = array_string;
    for (i = 0; i < new_length; ++i)
    {
        if (beginptr == NULL)
        {
            (*double_array)[i] = 0;
        }
        else
        {
            (*double_array)[i] = (double) strtof (beginptr, &endptr);
            beginptr = endptr;
        }
    }
}

void fclaw_options_destroy_array(void* array)
{
    FCLAW_FREE (array);
}
