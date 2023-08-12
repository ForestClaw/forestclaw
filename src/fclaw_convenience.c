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

#include "fclaw2d_convenience.h"
#include "fclaw3d_convenience.h"


void fclaw_domain_list_levels (fclaw_domain_t * domain, int log_priority)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_list_levels(domain,log_priority);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_list_levels(domain,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}

void fclaw_domain_list_neighbors (fclaw_domain_t * domain, int log_priority)
{
    if(domain->dim == 2)
    {
        fclaw2d_domain_list_neighbors(domain,log_priority);
    }
    else if(domain->dim == 3)
    {
        fclaw3d_domain_list_neighbors(domain,log_priority);
    }
    else
    {
        SC_ABORT_NOT_REACHED ();
    }
}