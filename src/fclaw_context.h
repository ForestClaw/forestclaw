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

/**
 * @file
 * 
 * @brief This file contains the declaration of the fclaw_context_t structure and related functions.
 * 
 * The fclaw_context_t structure represents a context for functions that can be paused and resumed.
 * It provides functions to retrieve and save values to the context object.
 * 
 * Example usage:
 * 
 * // Get the context object from the global context
 * fclaw_context_t *context = fclaw_context_get(glob, "my_context");
 * 
 * // Get an integer value from the context object
 * int value;
 * fclaw_context_get_int(context, "my_value", &value, 10);
 * 
 * // Get a double value from the context object
 * double value;
 * fclaw_context_get_double(context, "my_value", &value, 3.14);
 * 
 * // Save values to the context object
 * // This should be called right before exit points in a function
 * fclaw_context_save(context);
 *
 * This is designed with some restrictions in mind:
 * - The context object shoulde be retrieved at the beginning of a function with fclaw_context_get.
 * - The values should be retrieved with fclaw_context_get_int or fclaw_context_get_double.
 *    - These should be called in the same manner every time the resumable function is called.
 * - The values should be saved with fclaw_context_save right before exit points in the function.
 *    - the save call will save all the variables pointed to in the get_double/int calls.
 * 
 */

#ifndef FCLAW_CONTEXT_H
#define FCLAW_CONTEXT_H

#include <fclaw_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/**
 * @brief Context structure
 */
typedef struct fclaw_context fclaw_context_t;

/**
 * @brief Get the context object from glob. If the context object does not exist, it is created.
 * 
 * If the context exists, calling the function will reset all pointers associated with the values.
 *
 * fclaw_context_save needs to be called before the context object can be retrieved again.
 * 
 * @param glob the global context
 * @return fclaw_context_t* the context object
 */
fclaw_context_t* fclaw_context_get(fclaw_global_t *glob, const char *name);

/**
 * @brief Get an integer value from the context object.
 *
 * This function retrieves an integer value from the context object. If the value does not exist in the context,
 * the current value remains unchanged.
 * 
 * @param context the context object
 * @param name the name of the value
 * @param value a pointer to the value
 */
void fclaw_context_get_int(fclaw_context_t *context, 
                           const char *name,
                           int *value);

/**
 * @brief Get a double value from the context object. 
 *
 * This function retrieves an integer value from the context object. If the value does not exist in the context,
 * the current value remains unchanged.
 *
 * @param context the context object
 * @param name the name of the value
 * @param value a pointer to the value
 */
void fclaw_context_get_double(fclaw_context_t *context, 
                              const char *name, 
                              double *value);

/**
 * @brief Save values to the context object. Should be called right before exit points in a function.
 * 
 * This will get the values from the pointers provided in the get functions. 
 *
 * For values that were not retrieved between save calls, the old existing value will be kept.
 *
 * @param context the context object to save to
 */
void fclaw_context_save(fclaw_context_t *context);

/**
 * @brief Initializes the vtable needed for packing context objects
 * 
 * @param glob the global context
 */
void fclaw_context_vtable_initialize(fclaw_global_t *glob);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !FCLAW_MATH_H */
