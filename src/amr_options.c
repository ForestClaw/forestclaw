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

#include "amr_options.h"

void
amr_options_register (sc_options_t * opt, amr_options_t * amropt)
{
    sc_options_add_int (opt, 0, "mx-leaf", &amropt->mx_leaf, 32,
                        "Subdivision of each patch in x");
    sc_options_add_int (opt, 0, "my-leaf", &amropt->my_leaf, 32,
                        "Subdivision of each patch in y");
    sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 4.,
                           "Final time");
    sc_options_add_string (opt, 0, "subcycling", &amropt->subcycling,
                           "T", "Subcyclingy type");
    sc_options_add_inifile (opt, 'F', "inifile",
                            "Read options from this file");
}
  
void
amr_options_parse (sc_options_t * opt, int argc, char ** argv, int log_priority)
{
  int                   retval;

  retval = sc_options_parse (sc_package_id, SC_LP_ERROR, opt, argc, argv);
  if (retval < 0) {
    sc_options_print_usage (sc_package_id, log_priority, opt, NULL);
    sc_abort_collective ("Option parsing failed");
  }
  sc_options_print_summary (sc_package_id, log_priority, opt);
  if (sc_is_root ()) {
    retval = sc_options_save (sc_package_id, SC_LP_ERROR, opt,
                              "clawez.data.used");
    SC_CHECK_ABORT (!retval, "Option save failed");
  }
}
