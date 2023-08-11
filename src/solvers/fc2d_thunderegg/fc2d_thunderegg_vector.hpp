/*
  Copyright (c) 2019-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton, Grady Wright
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
 * Routines for interfaceing with ThunderEgg vectors
 */

#include <ThunderEgg/Vector.h>


/* Avoid circular dependencies */
struct fclaw_patch;
struct fclaw2d_domain;
struct fclaw2d_global;

/**
 * @brief Choice for patch data
 */
typedef enum fc2d_thunderegg_data_choice
{
    /** @brief RHS patch data */
    RHS=0,
    /** @brief soln patch data */
    SOLN,
    /** @brief soln patch data */
    STORE_STATE,
}  fc2d_thunderegg_data_choice_t;

/**
 * @brief Get a thunderegg vector that is a view to forestclaw data
 * 
 * @param glob the global context
 * @param data_choice the data choice
 * @return ThunderEgg::Vector<2> the vector
 */
ThunderEgg::Vector<2> fc2d_thunderegg_get_vector(struct fclaw2d_global *glob, fc2d_thunderegg_data_choice_t data_choice);

/**
 * @brief Get a thunderegg vector that is a view to forestclaw data
 * 
 * @param glob the global context
 * @param data_choice the data choice
 * @return ThunderEgg::Vector<2> the vector
 */
void fc2d_thunderegg_store_vector(struct fclaw2d_global *glob, fc2d_thunderegg_data_choice_t data_choice, const ThunderEgg::Vector<2>& vec);

