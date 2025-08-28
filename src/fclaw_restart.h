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

#ifndef FCLAW_RESTART_H
#define FCLAW_RESTART_H

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw_global;  /* This is a hack !! */


/**
 * @brief Restarts the forestclaw simulation from a restart file.
 *
 * This function is responsible for restarting the simulation from a previously saved state.
 * This is meant to be called from fclaw_initalize()
 *
 * This will print an error message and exit for any error.
 *
 * @param glob The global context
 * @param restart_filename The filename of the restart file.
 * @param partition_filename The filename of the partition file.
 */
void
fclaw_restart_from_file (struct fclaw_global * glob,
                         const char* restart_filename,
                         const char* partition_filename);

void
fclaw_reinitialize_from_file (struct fclaw_global * glob,
                              const char* restart_filename,
                              const char* partition_filename);

/**
 * @brief Tests the restart functionality by reading data from a restart file.
 *
 * This function is used to do an "in memory" restart test.
 * It deletes the current domain and creates a new domain from the restart file,
 * and overwrites the values in glob with the values from the restart file.
 *
 * @param glob The pointer to the global context
 * @param restart_filename The filename of the restart file to read from.
 * @param partition_filename The filename of the partition file to read from.
 */
void
fclaw_restart_test_from_file (struct fclaw_global * glob,
                              const char* restart_filename,
                              const char* partition_filename);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
