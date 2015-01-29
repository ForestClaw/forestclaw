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
 * It comes with its own option section "Blackbox"
 * containing the integer option "factor".
 */
typedef struct dummy_blackbox dummy_blackbox_t;

/** Create a new blackbox multiplicator object and register its options.
 * It will be deallocated when \ref fclaw_app_destroy is called.
 * \param [in] a                Allocated application that will deallocate
 *                              this object when it is destroyed.
 * \param [in] factor           The initial value for the multiplier.
 *                              We only consider multipliers in 0 and 10
 *                              inclusive valid values, even though this
 *                              function does not check it.  It is only
 *                              checked during option parsing.
 * \return                      Valid blackbox object.
 */
dummy_blackbox_t *dummy_blackbox_new_register (fclaw_app_t * a, int factor);

/** Perform the main operation of this object, integer multiplication.
 * \param [in] bbox             The blackbox multiplicator object.
 * \param [in] input            Arbitrary integer value.
 * \return                      The \b input multiplied with this box's factor.
 */
int dummy_blackbox_multiply (dummy_blackbox_t * bbox, int input);
