/*
  Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "swirl_user.h"
#include <fclaw2d_file.h>

#include "../all/advection_user.h"
#include <p4est_wrap.h>         /* just temporary for testing */
#include <fclaw2d_convenience.h>

#define FCLAW_SWIRL_IO_DEMO 0

static
void create_domain(fclaw_global_t *glob)
{
    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    /* Mapped, multi-block domain */
    fclaw_domain_t *domain = 
          fclaw_domain_new_unitsquare(glob->mpicomm, 
                                        fclaw_opt->minlevel);
    /* Create "empty" mapping */
    fclaw_map_context_t* cont = fclaw_map_new_nomap();

    /* Store domain in the glob */
    fclaw_global_store_domain(glob, domain);

    /* Map unit square to disk using mapc2m_disk.f */
    fclaw_map_store (glob, cont);

    /* Print out some info */
    fclaw_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);
}

#if FCLAW_SWIRL_IO_DEMO
static void
check_fclaw2d_file_error_code (int errcode, const char *str)
{
    int reslen, retval;
    char err_str[sc_MPI_MAX_ERROR_STRING];

    if (errcode != FCLAW2D_FILE_ERR_SUCCESS)
    {
        /* In case of a not successful fclaw2d_file function call we always
         * close the file if applicable and deallocate the file context.
         */
        /* examine the error code */
        retval = fclaw2d_file_error_string (errcode, err_str, &reslen);
        /* check for error in the error string function */
        SC_CHECK_ABORTF (!retval, "%s: error string function not successful",
                         str);
        SC_ABORTF ("%s: %*.*s", str, reslen, reslen, err_str);
    }
}
#endif

static
void run_program(fclaw_global_t* glob)
{
#if FCLAW_SWIRL_IO_DEMO
    int i;
    int errcode, retval;
    fclaw2d_file_context_t *fc;
    char read_user_string[FCLAW2D_FILE_USER_STRING_BYTES + 1];
    sc_array_t block_arr, field_arr, read_arr, *current_arr;
    int64_t test_int = 12;
    char *data, *local_arr_data;
    fclaw2d_domain_t *read_domain;
#endif

    /* Initialize virtual table for ForestClaw */
    fclaw_vtables_initialize(glob);

    /* Initialize virtual tables for solvers */
    const user_options_t *user_opt = swirl_get_options(glob);
    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_solver_initialize(glob);
    }
    else if (user_opt->claw_version == 5)
    {
        fc2d_clawpack5_solver_initialize(glob);
    }

    swirl_link_solvers(glob);

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fclaw_initialize(glob);
    fclaw_run(glob);

#if FCLAW_SWIRL_IO_DEMO
    /* Example usage of forestclaw file functions. This is just for
     * demonstration purposes. For an actual restart functionality in ForestClaw
     * the workflow must be extended by providing buffers with the required
     * data and the functions may be called at a more suitable place.
     */
    /* create a file, which is open for further writing */
    /* the passed domain is written to the file */
    fc = fclaw2d_file_open_write ("swirl_io_test", "ForestClaw data file",
                                  glob->domain->d2, &errcode);
    check_fclaw2d_file_error_code (errcode, "file open write");

#if 0
    retval = fclaw2d_file_write_partition ("swirl_io_test_partition",
                                           "Test partition write",
                                           glob->domain->d2, &errcode);
    check_fclaw2d_file_error_code (errcode, "file write partition");
#endif

    /* write a block to the file */
    /* Initialize a sc_array with one element and the element size equals
     * to the number of bytes of the block section.  */
    sc_array_init_data (&block_arr, &test_int, sizeof (int64_t), 1);
    fc = fclaw2d_file_write_block (fc, "Test block", block_arr.elem_size,
                                   &block_arr, &errcode);
    check_fclaw2d_file_error_code (errcode, "file write block");

    /* write an array associated to the domain to the file */
    /* we write non-contiguous data to demonstrate how to assemble the array */
    /* To this end, we initialize a sc_array with the number of local patches
     * of the domain passed to \ref fclaw2d_file_open_write and sizeof (sc_array_t)
     * as element size.
     */
    sc_array_init_size (&field_arr, sizeof (sc_array_t),
                        glob->domain->local_num_patches);

    for (i = 0; i < glob->domain->local_num_patches; ++i)
    {
        /* Each sc_array in field_array represents one data entity
         * associated to a patch. That means we write the data
         * associated to the i-th local patch below.
         */
        current_arr = (sc_array_t *) sc_array_index (&field_arr, i);
        /* To not allocate data but point to already allocated data
         * use \ref sc_array_init_data.
         */
        sc_array_init_size (current_arr, 3 * sizeof (char), 1);
        data = (char *) sc_array_index (current_arr, 0);
        data[0] = 'a';
        data[1] = 'b';
        data[2] = 'c';
    }

    fc = fclaw2d_file_write_array (fc, "Test array", 3 * sizeof (char),
                                   &field_arr, &errcode);
    check_fclaw2d_file_error_code (errcode, "file write array");
    /* free the local array data */
    for (i = 0; i < glob->domain->local_num_patches; ++i)
    {
        current_arr = (sc_array_t *) sc_array_index (&field_arr, i);
        sc_array_reset (current_arr);
    }
    sc_array_reset (&field_arr);
    retval = fclaw2d_file_close (fc, &errcode);
    check_fclaw2d_file_error_code (errcode, "file close 1");
    FCLAW_EXECUTE_ASSERT_FALSE (retval);

    /* open the file for reading */
    /* the domain stored in the file is read to read_domain */
    fc = fclaw2d_file_open_read ("swirl_io_test", read_user_string,
                                 glob->domain->mpicomm, NULL, &read_domain,
                                 &errcode);
    check_fclaw2d_file_error_code (errcode, "file open read");
    fclaw_global_productionf ("Opened file with user string: %s\n",
                              read_user_string);

    /* read a block from the file */
    test_int = -1;
    fc = fclaw2d_file_read_block (fc, read_user_string, sizeof (int64_t),
                                  &block_arr, &errcode);
    check_fclaw2d_file_error_code (errcode, "file read block");
    fclaw_global_productionf ("Read block with user string: %s\n",
                              read_user_string);
    FCLAW_ASSERT (test_int == 12);

    /* read an array from the file */
    /* For reading array data we need to pass a sc_array with element
     * equals to sizeof (sc_array_t). The sc_array will be resized by
     * \ref fclaw2d_file_read_array. Each entry of the output sc_array
     * is an sc_array with one element and element size equals to the
     * patch data size.
     */
    sc_array_init (&read_arr, sizeof (sc_array_t));
    fc = fclaw2d_file_read_array (fc, read_user_string, 3 * sizeof (char),
                                  &read_arr, &errcode);
    check_fclaw2d_file_error_code (errcode, "file write array");
    fclaw_global_productionf ("Read array with user string: %s\n",
                              read_user_string);
    /* check read array */
    for (i = 0; i < read_domain->local_num_patches; ++i)
    {
        current_arr = (sc_array_t *) sc_array_index_int (&read_arr, i);
        /* local_arr_data is the pointer to the actual data of the
         * i-th array entry
         */
        local_arr_data = (char *) sc_array_index (current_arr, 0);
        /* check array entry */
        FCLAW_ASSERT (local_arr_data[0] == 'a');
        FCLAW_ASSERT (local_arr_data[1] == 'b');
        FCLAW_ASSERT (local_arr_data[2] == 'c');
    }
    for (i = 0; i < read_domain->local_num_patches; ++i)
    {
        current_arr = (sc_array_t *) sc_array_index (&read_arr, i);
        sc_array_reset (current_arr);
    }
    sc_array_reset (&read_arr);

    /* sanity check of read domain */
    FCLAW_ASSERT (p4est_checksum (((p4est_wrap_t *) read_domain->pp)->p4est)
                  ==
                  p4est_checksum (((p4est_wrap_t *) glob->domain->d2->pp)->
                                  p4est));

    fclaw2d_domain_destroy (read_domain);

    retval = fclaw2d_file_close (fc, &errcode);
    check_fclaw2d_file_error_code (errcode, "file close 2");
    FCLAW_EXECUTE_ASSERT_FALSE (retval);
#endif

    fclaw_finalize(glob);
}

int
main (int argc, char **argv)
{
    /* Initialize application */
    fclaw_app_t *app = fclaw_app_new (&argc, &argv, NULL);

    /* Options */
    user_options_t              *user_opt;
    fclaw_options_t             *fclaw_opt;
    fclaw2d_clawpatch_options_t *clawpatch_opt;
    fc2d_clawpack46_options_t   *claw46_opt;
    fc2d_clawpack5_options_t    *claw5_opt;

    /* Create new options packages */
    fclaw_opt =                   fclaw_options_register(app,  NULL,        "fclaw_options.ini");
    clawpatch_opt =   fclaw2d_clawpatch_options_register(app, "clawpatch",  "fclaw_options.ini");
    claw46_opt =        fc2d_clawpack46_options_register(app, "clawpack46", "fclaw_options.ini");
    claw5_opt =          fc2d_clawpack5_options_register(app, "clawpack5",  "fclaw_options.ini");
    user_opt =                    swirl_options_register(app,               "fclaw_options.ini");

    /* Read configuration file(s) and command line, and process options */
    int first_arg;
    fclaw_exit_type_t vexit;
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    /* Run the program */
    if (!vexit)
    {
        /* Create global structure which stores the domain, timers, etc */
        int size, rank;
        sc_MPI_Comm mpicomm = fclaw_app_get_mpi_size_rank (app, &size, &rank);
        fclaw_global_t *glob = fclaw_global_new_comm (mpicomm, size, rank);

        /* Store option packages in glob */
        fclaw_options_store           (glob, fclaw_opt);
        fclaw2d_clawpatch_options_store (glob, clawpatch_opt);
        fc2d_clawpack46_options_store   (glob, claw46_opt);
        fc2d_clawpack5_options_store    (glob, claw5_opt);
        swirl_options_store             (glob, user_opt);

        //char buffer[fclaw2d_global_packsize(glob)];
        //fclaw2d_global_pack(glob,buffer);
        //fclaw2d_global_t* glob2;
        //fclaw2d_global_unpack(buffer, &glob2);

        /* Create and store domain */
        create_domain(glob);

        run_program(glob);

        fclaw_global_destroy(glob);
        //fclaw2d_global_destroy(glob2);
    }

    fclaw_app_destroy (app);

    return 0;
}
