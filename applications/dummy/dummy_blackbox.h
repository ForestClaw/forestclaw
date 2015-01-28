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

/** \file dummy_blackbox.h
 * This file is a dummy package to do integer multiplication.
 * Its purpose is to show how a package can hide their options from main.
 */

#include <fclaw_base.h>

/** The blackbox multiplicator is an paque object.
 * Right now it is setup to call its constructor, register its options,
 * and eventually destroy the object.  Care must be taken to destroy
 * this object after calling \ref fclaw_app_destroy.
 *
 * We might also change this interface to hide the new and destroy functions
 * altogether and handle them through the options virtual table.
 * This would shorten the main program by two lines and make the workings
 * of this blackbox object even more obscure.  On the plus side, above
 * requirements on the order of \ref dummy_blackbox_destroy would vanish.
 */
typedef struct dummy_blackbox dummy_blackbox_t;

/** This virtual table can be used in registering options.
 * This will add an integer parameter -f / --factor to selecet a multiplier.
 * We only consider multipliers in 0 and 10 inclusive valid values.
 */
extern const fclaw_app_options_vtable_t *dummy_blackbox_vt;

/** Create a now blackbox multiplicator object.
 * \param [in] factor           The initial value for the multiplier.
 *                              We only consider multipliers in 0 and 10
 *                              inclusive valid values.
 * \return                      Valid blackbox object.
 */
dummy_blackbox_t *dummy_blackbox_new (int factor);

/** Destroy a blackbox object.
 * \note It is important to call this after \ref fclaw_app_destroy,
 *       since that function activates the options destructors that
 *       might access this blackbox object.
 * \param [in,out] bbox         The object is deallocated.
 */
void dummy_blackbox_destroy (dummy_blackbox_t * bbox);

/** Perform the main operation of this object, integer multiplication.
 * \param [in] bbox             The blackbox multiplicator object.
 * \param [in] input            Arbitrary integer value.
 * \return                      The \b input multiplied with this box's factor.
 */
int dummy_blackbox_multiply (dummy_blackbox_t * bbox, int input);
