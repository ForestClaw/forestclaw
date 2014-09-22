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

#include "amr_forestclaw.H"
#include "amr_waveprop.H"
#include "test_parms_user.H"

#if 0
#define TEST_ALLOC_ZERO SC_ALLOC_ZERO
#define TEST_FREE       SC_FREE
#else
#define TEST_ALLOC_ZERO FCLAW2D_ALLOC_ZERO
#define TEST_FREE       FCLAW2D_FREE
#endif

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif


amr_waveprop_parms_t*  test_parms_new(sc_options_t *opt)
{
    /* This is a good example of a place where a [Section] would be nice.
       Something like [waveprop].  See fclaw_defaults.ini */

    amr_waveprop_parms_t *waveprop_parms;

    waveprop_parms = TEST_ALLOC_ZERO(amr_waveprop_parms_t, 1);

    /* Array of SpaceDim many values, with no defaults is set to all 0's */
    amr_options_add_int_array (opt, 0, "order", &waveprop_parms->order_string, NULL,
                               &waveprop_parms->order, SpaceDim,
                               "Normal and transverse orders");

   sc_options_load (sc_package_id, SC_LP_ALWAYS, opt, "fclaw2d_waveprop.ini");

   // amr_waveprop_postprocess_parms(waveprop_parms);

   return waveprop_parms;

}


void test_parms_postprocess_parms(amr_waveprop_parms_t* waveprop_parms)
{
    /* -----------------------------------------------------------------------
       Some post-processing of arrays
       ------------------------------------------------------------------------ */
    amr_options_convert_int_array (waveprop_parms->order_string, &waveprop_parms->order,
                                   SpaceDim);


}


void test_parms_delete(amr_waveprop_parms_t* waveprop_parms)
{
    SC_FREE(waveprop_parms->order);
    TEST_FREE(waveprop_parms);
    waveprop_parms = NULL;
}


#ifdef __cplusplus
#if 0
{
#endif
}
#endif
