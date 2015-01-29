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

#include "dummy_blackbox.h"

/** This is the private declaration of the blackbox object. */
struct dummy_blackbox
{
    int factor;                 /**< A factor used for multiplication. */
};

static void *
dummy_blackbox_register (fclaw_app_t * a, void *package,
                         sc_options_t * options)
{
    dummy_blackbox_t *bbox = (dummy_blackbox_t *) package;

    FCLAW_ASSERT (bbox != NULL);

    sc_options_add_int (options, 'f', "factor", &bbox->factor,
                        bbox->factor, "Integer multiplier");

    return NULL;
}

static fclaw_exit_type_t
dummy_blackbox_check (fclaw_app_t * a, void *package, void *registered)
{
    dummy_blackbox_t *bbox = (dummy_blackbox_t *) package;

    FCLAW_ASSERT (bbox != NULL);
    FCLAW_ASSERT (registered == NULL);

    if (bbox->factor < 0 || bbox->factor > 10)
    {
        fclaw_global_errorf
            ("The value for the blackbox factor must be between 0 and 10\n");
        return FCLAW_EXIT_ERROR;
    }

    return FCLAW_NOEXIT;
}

static void
dummy_blackbox_destroy (fclaw_app_t * a, void *package, void *registered)
{
    dummy_blackbox_t *bbox = (dummy_blackbox_t *) package;

    FCLAW_ASSERT (bbox != NULL);
    FCLAW_ASSERT (registered == NULL);

    FCLAW_FREE (bbox);
}

static const fclaw_app_options_vtable_t dummy_blackbox_vt = {
    dummy_blackbox_register,
    NULL,
    dummy_blackbox_check,
    dummy_blackbox_destroy
};

dummy_blackbox_t *
dummy_blackbox_new_register (fclaw_app_t * a, int factor)
{
    dummy_blackbox_t *bbox = FCLAW_ALLOC (dummy_blackbox_t, 1);

    bbox->factor = factor;
    fclaw_app_options_register (a, "Blackbox", NULL, &dummy_blackbox_vt,
                                bbox);

    return bbox;
}

int
dummy_blackbox_multiply (dummy_blackbox_t * bbox, int input)
{
    FCLAW_ASSERT (bbox != NULL);

    return bbox->factor * input;
}
