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

#include <sc_statistics.h>
#include <fclaw2d_global.h>

#include <fclaw2d_diagnostics.h>
#include <fclaw2d_options.h>
#include <fclaw2d_map.h>
#include <fclaw2d_domain.h>
#include <fclaw2d_forestclaw.h>


static void write_output(FILE* file, const char* format, ...)
{

    va_list args;
    va_start(args, format);
    int fprintf_value = vfprintf(file, format, args);
    va_end(args);

    if(fprintf_value < 0)
    {
        fclaw_global_essentialf("Failed to write to output file\n");
        exit(1);
    }
}

static int read_input_and_check(FILE* file, const char* name, int actual)
{
    int name_length = strlen(name);
    char* format = FCLAW_ALLOC(char, name_length+4+1);
    strncpy(format, name, name_length + 1);
    strcat(format, ",%d\n");

    int expected;
    int fprintf_value = fscanf(file, format, &expected);

    if(fprintf_value < 0)
    {
        fclaw_global_essentialf("Failed to write to output file\n");
        exit(1);
    }

    int failure = 0;
    if(expected != actual)
    {
        fclaw_global_essentialf("Expected %s to be %d, but was %d\n", name, expected, actual);
        failure = 1;
    }

    FCLAW_FREE(format);

    return failure;
}

/**
 * @brief Output expected values for regression testing in a fortran style output
 * 
 * @param glob the global structure
 * @param filename the base name of the output file
 */
static void 
output_expected_values(fclaw2d_global_t* glob, const char* filename)
{
    if(glob->mpirank == 0)
    {
        // append .expected to the filename
        char* expected_filename = FCLAW_ALLOC(char, strlen(filename) + 10);
        strcpy(expected_filename, filename);
        strcat(expected_filename, ".expected");

        // open the file
        FILE* file = fopen(expected_filename, "w");
        if (file == NULL)
        {
            fclaw_global_essentialf("Failed to open expected values file %s\n", expected_filename);
            exit(FCLAW_EXIT_ERROR);
        }

        //write amr_dvance in fortran style output
        write_output(file, "count_amr_advance,%d\n", glob->count_amr_advance);
        write_output(file, "global_num_patches,%d\n", glob->domain->global_num_patches);
        write_output(file, "count_amr_new_domain,%d\n", glob->count_amr_new_domain);

        //close file
        fclose(file);

        FCLAW_FREE(expected_filename);

        // open file to get expected results
        file = fopen(filename, "r");
        if (file == NULL)
        {
            fclaw_global_essentialf("Failed to open regressions file %s\n", expected_filename);
            exit(FCLAW_EXIT_ERROR);
        }

        int num_failures = 0;
        num_failures += read_input_and_check(file, "count_amr_advance", glob->count_amr_advance);
        num_failures += read_input_and_check(file, "global_num_patches", glob->domain->global_num_patches);
        num_failures += read_input_and_check(file, "count_amr_new_domain", glob->count_amr_new_domain);

        if(num_failures > 0){
            exit(1);
        }

        // close file
        fclose(file);
    }
}
/* ------------------------------------------------------------------
   Public interface
   ---------------------------------------------------------------- */

void fclaw2d_finalize(fclaw2d_global_t* glob)
{
    const fclaw_options_t *gparms = fclaw2d_get_options(glob);

    fclaw_global_essentialf("Finalizing run\n");
    fclaw2d_diagnostics_finalize(glob);
    if (glob->cont != NULL) {
        fclaw2d_map_destroy(glob->cont);
    }
    fclaw2d_domain_barrier (glob->domain);

    if (gparms->report_timing)
    {
        if (gparms->outstyle > 0)
        {
            /* Only call this if we have taken time steps.  For time-independent problems, we
               probably need a different report (no "amr_advance_steps") */
            fclaw2d_timer_report(glob);
        }
        else
        {
            fclaw_global_essentialf("Timing reports not generated for outstyle=0\n");
        }
    }
    if (gparms->regression_check)
    {
        output_expected_values(glob, gparms->regression_check);
    }
    fclaw2d_domain_reset(glob);
}
