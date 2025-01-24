/*
Copyright (c) 2012-2024 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw_restart.h>
#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_vtable.h>
#include <fclaw_patch.h>
#include <fclaw_domain.h>
#include <fclaw_exchange.h>
#include <fclaw_file.h>
#include <fclaw_convenience.h>
#include <fclaw_regrid.h>
#include <fclaw_partition.h>
#include <fclaw_forestclaw.h>
#include <fclaw_output.h>
#include <iniparser.h>

/**
 * @brief Check the error code and print an error message if not successful
 */
#define CHECK_ERROR_CODE(refine_dim, errcode, str) \
do { \
    int reslen, retval; \
    char err_str[sc_MPI_MAX_ERROR_STRING]; \
    if (errcode != FCLAW_FILE_ERR_SUCCESS) \
    { \
        retval = fclaw_file_error_string (refine_dim, errcode, err_str, &reslen); \
        if(retval != 0) \
        { \
            fclaw_global_errorf ("%s: error string function not successful", str); \
        } \
        fclaw_global_errorf("%s: %s\n", str, err_str); \
    } \
} while(0)

/**
 * @brief Check the error code and abort with message if not successful
 */
#define CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, str) \
do { \
    int reslen, retval; \
    char err_str[sc_MPI_MAX_ERROR_STRING]; \
    if (errcode != FCLAW_FILE_ERR_SUCCESS) \
    { \
        retval = fclaw_file_error_string (refine_dim, errcode, err_str, &reslen); \
        SC_CHECK_ABORTF (!retval, "%s: error string function not successful", str); \
        SC_ABORTF ("%s: %*.*s", str, reslen, reslen, err_str); \
    } \
} while(0)


/**
 * @brief Get the ini file for the options as a string
 * 
 * @param glob the global context
 * @return char* the ini file as a string, NULL if an error occurs. Needs to be freed by the caller
 */
static char*
get_used_ini(fclaw_global_t * glob)
{
    char *buffer = NULL;
    long length = 0;
    if(glob->mpirank == 0)
    {
        sc_options_t * options = (sc_options_t*) fclaw_global_get_attribute(glob, "fclaw_options");
        if(options == NULL)
        {
            fclaw_abortf("fclaw_restart.c: Cannot find sc options structure in glob."
                         "                 Make sure the fclaw_global_new(app) constructor is being used in main.\n");
        }

        int retval = sc_options_save (fclaw_get_package_id (),
                                      FCLAW_VERBOSITY_ERROR, 
                                      options, 
                                      "fclaw_options.ini.used");
        if(retval != 0)
        {
            fclaw_global_productionf("fclaw_restart.c: Error saving options\n");
            return NULL;
        }

        // read the entire file into a string
        FILE *file = fopen("fclaw_options.ini.used", "r");
        if (file == NULL)
        {
            printf("fclaw_restart.c: Cannot open file\n");
            return NULL;
        }

        fseek(file, 0, SEEK_END);
        length = ftell(file);
        fseek(file, 0, SEEK_SET);

        buffer = FCLAW_ALLOC(char, length + 1);
        if (buffer)
        {
            fread(buffer, 1, length, file);
        }
        fclose(file);

        buffer[length] = '\0';
    }
    //broadcast the lenth
    sc_MPI_Bcast(&length, 1, sc_MPI_LONG, 0, glob->mpicomm);
    //allocate the buffer on other ranks
    if(glob->mpirank != 0)
    {
        buffer = FCLAW_ALLOC(char, length + 1);
    }
    FCLAW_ASSERT(buffer != NULL);
    //broadcast the string
    sc_MPI_Bcast(buffer, length+1, sc_MPI_CHAR, 0, glob->mpicomm);

    return buffer;
    

}

/**
 * @brief Check if option should be skipped for comparison
 * 
 * @param fclaw_opt_sections the sections containing the core options
 * @param section the section name
 * @param key the key name
 * @return int true if the option should be skipped, false otherwise
 */
static int
skip_option(sc_keyvalue_t *fclaw_opt_sections, const char* section, const char* key)
{
    int skip = 0;
    if(sc_keyvalue_get_int(fclaw_opt_sections, section, 0))
    {
        //skip the section: part of the key
        const char* key_substr = key + strlen(section) + 1;
        if(strcmp(key_substr, "restart-file") == 0)
        {
            skip = 1;
        }
        else if(strcmp(key_substr, "partition-file") == 0)
        {
            skip = 1;
        }
        else if(strcmp(key_substr, "restart") == 0)
        {
            skip = 1;
        }
        else if(strcmp(key_substr, "checkpoint") == 0)
        {
            skip = 1;
        }
    }
    return skip;
}

/**
 * @brief compare two ini dictionaries
 * 
 * @param fclaw_opt_secitons the sections containing the core options
 * @param expected the expected dictionary
 * @param actual the actual dictionary
 * @return int the number of differences found
 */
static int
compare_dictionaries(sc_keyvalue_t *fclaw_opt_secitons,
                     dictionary *expected, 
                     dictionary *actual)
{
    int num_differences = 0;
    int nsec = iniparser_getnsec(expected);
    for (int i_sec = 0; i_sec < nsec; i_sec++)
    {
        char* section = iniparser_getsecname(expected, i_sec);
        if(strcmp(section, "arguments") != 0)
        {
            if(iniparser_find_entry(actual,section))
            {
                int nkey = iniparser_getsecnkeys(expected, section);
                char** keys = iniparser_getseckeys(expected, section);

                for (int i_key = 0; i_key < nkey; i_key++)
                {
                    char* key = keys[i_key];
                    if(skip_option(fclaw_opt_secitons, section, key))
                    {
                        continue;
                    }

                    if(iniparser_find_entry(actual, key))
                    {
                        char* expected_value = iniparser_getstring(expected, key, NULL);
                        char* actual_value = iniparser_getstring(actual, key, NULL);
                        if(strcmp(expected_value, actual_value) != 0)
                        {
                            fclaw_global_productionf("restarting with mismatched value for [%s]. Value in checkpoint [%s], value being run with [%s].\n", key, expected_value, actual_value);
                            num_differences++;
                        }

                    }
                    else
                    {
                        fclaw_global_productionf("Checkpoint file has unused option %s.\n", keys[i_key]);
                        num_differences++;
                    }
                }
                free (keys);
            }
            else
            {
                fclaw_global_productionf("Checkpoint file has unused section [%s].\n", section);
                num_differences++;
            }
        }
    }
    // do inverse and look for options in actual that are not in expected
    nsec = iniparser_getnsec(actual);
    for(int i_sec = 0; i_sec < nsec; i_sec++)
    {
        char* section = iniparser_getsecname(actual, i_sec);
        if(strcmp(section, "arguments") != 0)
        {
            if(!iniparser_find_entry(expected, section))
            {
                fclaw_global_productionf("Checkpoint file is missing section [%s].\n", section);
                num_differences++;
            }
            else
            {
                int nkey = iniparser_getsecnkeys(actual, section);
                char** keys = iniparser_getseckeys(actual, section);

                for(int i_key = 0; i_key < nkey; i_key++)
                {
                    char* key = keys[i_key];
                    if(skip_option(fclaw_opt_secitons, section, key))
                    {
                        continue;
                    }

                    if(!iniparser_find_entry(expected, key))
                    {
                        fclaw_global_productionf("Checkpoint file is missing option %s.\n", keys[i_key]);
                        num_differences++;
                    }
                }
                free(keys);
            }
        }
    }
    return num_differences;
}

/**
 * @brief Compare provided options with the options stored in the checkpoint
 * 
 * @param glob the global context
 * @param checkpoint_ini the ini file from the checkpoint
 */
static void
check_options(fclaw_global_t * glob, const char* checkpoint_ini)
{
    if(glob->mpirank == 0)
    {
        fclaw_global_productionf("\n");
        fclaw_global_productionf("=========== Comparing with options stored in checkpoint ===========\n");
        fclaw_global_productionf("\n");

        //save the used ini file
        sc_options_t * options = 
            (sc_options_t*) fclaw_global_get_attribute(glob, "fclaw_options");
        int retval = sc_options_save (fclaw_get_package_id (),
                                      FCLAW_VERBOSITY_ERROR, 
                                      options, 
                                      "fclaw_options.ini.used");
        
        if(retval != 0)
        {
            fclaw_abortf("fclaw_restart.c: Error saving current options\n");
        }

        //save the restart ini file to fclaw_options.ini.checkpoint
        FILE *file = fopen("fclaw_options.ini.checkpoint", "w");
        if (file == NULL)
        {
            printf("Cannot open file\n");
            return;
        }
        fprintf(file, "%s", checkpoint_ini);
        fclose(file);

        fclaw_global_productionf("Checkpoint options have been saved to fclaw_options.ini.checkpoint\n");
        fclaw_global_productionf("\n");

        dictionary *actual = iniparser_load("fclaw_options.ini.used");
        dictionary *expected = iniparser_load("fclaw_options.ini.checkpoint");

        sc_keyvalue_t *fclaw_opt_sections 
            = (sc_keyvalue_t*) fclaw_global_get_attribute(glob, "fclaw_opt_sections");
        int num_differences = compare_dictionaries(fclaw_opt_sections, expected, actual);
        
        iniparser_freedict(actual);
        iniparser_freedict(expected);

        if(num_differences == 0)
        {
            fclaw_global_productionf("No differences found.\n");
        }

        fclaw_global_productionf("\n");
        fclaw_global_productionf("===================================================================\n");
        fclaw_global_productionf("\n");
    }
}

static void
free_used_ini(void* data)
{
    char* buffer = (char*) data;
    FCLAW_FREE(buffer);
}

/**
 * @brief Check that the user string matches an expected string.
 * The user string returned by fclaw_file functions is padded with spaces,
 * so we need to check that the strings match up to the first space and that 
 * the rest of the string is just spaces.
 * 
 * @param expected  the expected string
 * @param actual the actual string
 */
static
void check_user_string(const char* expected, const char* actual)
{
    if(strncmp(expected, actual, strlen(expected)) != 0)
    {
        fclaw_abortf("fclaw_restart.c: User string mismatch: %s != %s\n", expected, actual);
    }
    else
    {
        //also check that rest of the string is just spaces
        for(int i = strlen(expected); i < FCLAW_FILE_USER_STRING_BYTES; i++)
        {
            if(actual[i] != '\0')
            {
                break;
            }
            else if(actual[i] != ' ')
            {
                fclaw_abortf("fclaw_restart.c: User string mismatch: %s != %s\n", expected, actual);
            }
        }
    }
}

typedef struct pack_iter
{
    fclaw_global_t * glob;
    int curr_index;
    size_t size;
    sc_array_t* patches;
    int pointerno;
}pack_iter_t;


/**
 * @brief Set patch data from the checkpoint, user data is a pack_iter_t
 */
static void
set_patches(fclaw_domain_t * domain, fclaw_patch_t * patch, int blockno, int patchno,  void *user)
{
    pack_iter_t *user_data = (pack_iter_t*)user;
    sc_array_t *patches = user_data->patches;
    sc_array_t * current_arr = (sc_array_t *) sc_array_index (patches, user_data->curr_index);

    fclaw_build_mode_t build_mode = FCLAW_BUILD_FOR_UPDATE;

    if(user_data->pointerno == 0)
    {
	    fclaw_patch_build(user_data->glob, domain, patch, blockno, patchno,(void*) &build_mode);
    }

    void* data = fclaw_patch_checkpoint_get_pointer(user_data->glob, patch, blockno, patchno, user_data->pointerno);

    memcpy(data, sc_array_index(current_arr, 0), user_data->size);

    sc_array_reset(current_arr);

    user_data->curr_index++;
}

static
void restart (fclaw_global_t * glob,
              const char* restart_filename,
              const char* partition_filename,
              fclaw_timer_names_t timer)
{
    int refine_dim = glob->domain->refine_dim;
    fclaw_domain_reset(glob);

    int errcode;
    sc_array_t* partition = NULL;
    char user_string[FCLAW_FILE_USER_STRING_BYTES];
    if(partition_filename != NULL)
    {
        partition = sc_array_new(sizeof(p4est_gloidx_t));
        fclaw_file_read_partition(refine_dim, 
                                  partition_filename, 
                                  user_string, 
                                  glob->mpicomm, 
                                  partition, 
                                  &errcode);
        check_user_string("Partition", user_string);
        CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read_partition");
    }

    fclaw_file_context_t *fc 
        = fclaw_file_open_read (refine_dim,
                                restart_filename, 
                                user_string, 
                                glob->mpicomm, 
                                partition, 
                                &glob->domain, 
                                &errcode);
    check_user_string("ForestClaw checkpoint file", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart open_file");

    fclaw_domain_setup(glob, glob->domain);

    if(partition != NULL)
    {
        sc_array_destroy(partition);
    }

    sc_array_t array;

    //read the version major

    sc_array_init_size(&array, sizeof(int), 1);
    fc = fclaw_file_read_block(fc, user_string, sizeof(int), &array, &errcode);

    check_user_string("checkpoint_version_major", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read version_major");

    int version_major = *((int*) sc_array_index(&array, 0));

    sc_array_reset(&array);


    //read the version minor

    sc_array_init_size(&array, sizeof(int), 1);
    fc = fclaw_file_read_block(fc, user_string, sizeof(int), &array, &errcode);

    check_user_string("checkpoint_version_minor", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read version_minor");

    int version_minor = *((int*) sc_array_index(&array, 0));

    sc_array_reset(&array);

    //version check
    if(version_major != 1 || version_minor != 0)
    {
        fclaw_abortf("fclaw_restart.c: Incompatible checkpoint version: %d.%d\n", version_major, version_minor);
    }

    //read the length of the used_ini string

    sc_array_init_size(&array, sizeof(size_t), 1);
    fc = fclaw_file_read_block(fc, user_string, sizeof(size_t), &array, &errcode);


    check_user_string("used_ini_length", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read used_ini_length");

    size_t ini_length = *((size_t*) sc_array_index(&array, 0));
    sc_array_reset(&array);

    //read the used_ini string

    sc_array_init_size(&array, ini_length, 1);
    fc = fclaw_file_read_block(fc, user_string, ini_length, &array, &errcode);


    check_user_string("used_ini", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read used_ini");

    const char* used_ini = (const char*) sc_array_index(&array, 0);
    check_options(glob, used_ini);
    
    sc_array_reset(&array);

    //read the length of the glob buffer

    sc_array_init_size(&array, sizeof(size_t), 1);

    fc = fclaw_file_read_block(fc, user_string, sizeof(size_t), &array, &errcode);

    check_user_string("glob_size", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read globsize");

    size_t glob_packsize = *((size_t*) sc_array_index(&array, 0));

    sc_array_reset(&array);

    //read the glob buffer

    sc_array_init_size(&array, glob_packsize, 1);

    fc = fclaw_file_read_block(fc, user_string, glob_packsize, &array, &errcode);

    check_user_string("glob", user_string);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read glob buffer");

    fclaw_global_unpack((char *) sc_array_index(&array, 0), glob);

    sc_array_reset(&array);

    //read patch data

    int num_pointers = fclaw_patch_checkpoint_num_pointers(glob);
    size_t sizes[num_pointers];
    fclaw_patch_checkpoint_pointer_sizes(glob, sizes);
    const char* names[num_pointers];
    fclaw_patch_checkpoint_names(glob, names);

    for(int i = 0; i < num_pointers; i++)
    {
        sc_array_t *patches = sc_array_new_count(sizeof(sc_array_t), glob->domain->local_num_patches);
        pack_iter_t user;
        user.glob = glob;
        user.curr_index = 0;
        user.patches = patches;
        user.size = sizes[i];
        user.pointerno = i;

        fc = fclaw_file_read_array(fc, user_string, sizes[i], patches, &errcode);

        check_user_string(names[i], user_string);
        CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart read patches");

        fclaw_domain_iterate_patches(glob->domain, set_patches, &user);

        sc_array_destroy(patches);
    }

    fclaw_file_close(fc, &errcode);
    CHECK_ERROR_CODE_AND_ABORT(refine_dim, errcode, "restart close file");

    fclaw_initialize_domain_flags(glob);
    fclaw_exchange_setup(glob,timer);
    fclaw_regrid_set_neighbor_types(glob);
}

/**
 * @brief Get patch data for the checkpoint, user data is a pack_iter_t
 */
static void
get_patches(fclaw_domain_t * domain, fclaw_patch_t * patch, int blockno, int patchno,  void *user)
{
    pack_iter_t *user_data = (pack_iter_t*)user;
    sc_array_t *patches = user_data->patches;
    sc_array_t * current_arr = (sc_array_t *) sc_array_index (patches, user_data->curr_index);

    void* data = fclaw_patch_checkpoint_get_pointer(user_data->glob, patch, blockno, patchno, user_data->pointerno);

    sc_array_init_data(current_arr, data, user_data->size, 1);

    user_data->curr_index++;
}


static void
checkpoint_output_frame (fclaw_global_t * glob, int iframe)
{
    int refine_dim = glob->domain->refine_dim;

    char filename[BUFSIZ];
    char parition_filename[BUFSIZ];
    snprintf(filename, BUFSIZ, "fort_frame_%04d.checkpoint", iframe);
    snprintf(parition_filename, BUFSIZ, "fort_frame_%04d.partition", iframe);

    int errcode;
    fclaw_file_context_t *fc 
        = fclaw_file_open_write (filename, "ForestClaw checkpoint file",
                                 glob->domain, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "checkpoint open file");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    sc_array_t array;

    //dont know if we are going to bother with supporting older versions
    //of checkpoint files, but we will write the version number anyway

    //write version major

    int version_major = 1;
    sc_array_init_data(&array, &version_major, sizeof(int), 1);

    fc = fclaw_file_write_block(fc, "checkpoint_version_major", sizeof(int), &array, &errcode);
    CHECK_ERROR_CODE(refine_dim, errcode, "checkpoint write version_major");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    //write version minor

    int version_minor = 0;
    sc_array_init_data(&array, &version_minor, sizeof(int), 1);
    fc = fclaw_file_write_block(fc, "checkpoint_version_minor", sizeof(int), &array, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "checkpoint write version_minor");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    //write the used_ini length and string

    char* used_ini = 
        (char*) fclaw_global_get_attribute(glob, "fclaw_used_ini");
    if(used_ini == NULL)
    {
        used_ini = get_used_ini(glob);
        if(used_ini == NULL)
        {
            return;
        }
        fclaw_global_attribute_store(glob, 
                                     "fclaw_used_ini", 
                                     used_ini, 
                                     NULL, 
                                     free_used_ini);
    }

    size_t used_ini_length = strlen(used_ini)+1;
    sc_array_init_data(&array, &used_ini_length, sizeof(size_t), 1);

    fc = fclaw_file_write_block(fc, "used_ini_length", sizeof(size_t), &array, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "write used_ini_length");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    sc_array_init_data(&array, used_ini, used_ini_length, 1);
    fc = fclaw_file_write_block(fc, "used_ini", used_ini_length, &array, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "write used_ini");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    //write size of glob buffer

    size_t glob_packsize = fclaw_global_packsize(glob);
    sc_array_init_data(&array, &glob_packsize, sizeof(size_t), 1);

    fc = fclaw_file_write_block(fc, "glob_size", sizeof(size_t), &array, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "write globsize");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    //write glob buffer

    sc_array_t glob_buffer;
    sc_array_init_size(&glob_buffer, glob_packsize, 1);
    fclaw_global_pack(glob,(char *) sc_array_index(&glob_buffer, 0));

    fc = fclaw_file_write_block(fc, "glob", glob_packsize, &glob_buffer, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "write glob buffer");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }

    sc_array_reset(&glob_buffer);

    //write patch data

    int num_pointers = fclaw_patch_checkpoint_num_pointers(glob);
    size_t sizes[num_pointers];
    fclaw_patch_checkpoint_pointer_sizes(glob, sizes);
    const char* names[num_pointers];
    fclaw_patch_checkpoint_names(glob, names);
    for(int i = 0; i < num_pointers; i++)
    {
        sc_array_t *patches = sc_array_new_count(sizeof(sc_array_t), glob->domain->local_num_patches);
        pack_iter_t user;
        user.glob = glob;
        user.curr_index = 0;
        user.patches = patches;
        user.size = sizes[i];
        user.pointerno = i;
        fclaw_domain_iterate_patches(glob->domain, get_patches, &user);

    
        fc = fclaw_file_write_array(fc, names[i], sizes[i], patches, &errcode);
        CHECK_ERROR_CODE(refine_dim , errcode, "write patches");
        if(errcode != FCLAW_FILE_ERR_SUCCESS)
        {
            return;
        }

        for(int i = 0; i < glob->domain->local_num_patches; i++)
        {
            sc_array_t * current_arr = (sc_array_t *) sc_array_index (patches, i);
            sc_array_reset(current_arr);
        }
        sc_array_destroy(patches);
    }


    fclaw_file_close(fc, &errcode);
    fclaw_file_write_partition (parition_filename,
                                "Paritition",
                                glob->domain, &errcode);
    CHECK_ERROR_CODE(refine_dim , errcode, "close file");
    if(errcode != FCLAW_FILE_ERR_SUCCESS)
    {
        return;
    }
}

/* -----------------------------------------------------------------------
    Public interface
    -------------------------------------------------------------------- */

void
fclaw_restart_from_file (fclaw_global_t * glob,
                         const char* restart_filename,
                         const char* partition_filename)
{
    restart(glob, restart_filename, partition_filename, FCLAW_TIMER_INIT);
}

void fclaw_output_checkpoint(fclaw_global_t* glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);

    if(fclaw_opt->checkpoint)
    {
        fclaw_timer_start (&glob->timers[FCLAW_TIMER_OUTPUT]);

        checkpoint_output_frame(glob,iframe);

        fclaw_timer_stop (&glob->timers[FCLAW_TIMER_OUTPUT]);
    }
}