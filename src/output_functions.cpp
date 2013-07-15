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

#include "amr_utils.H"
#include "amr_output.H"

fclaw2d_output_functions_t* get_output_functions(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t* ddata = get_domain_data(domain);
    return ddata->output_functions;
}

void link_output_functions(fclaw2d_domain_t* domain,
                           fclaw2d_patch_write_header_t patch_write_header,
                           fclaw2d_patch_write_output_t patch_write_output)
{
    fclaw2d_output_functions_t *of = get_output_functions(domain);

    of->f_patch_write_header = patch_write_header;
    of->f_patch_write_output = patch_write_output;
}

void initialize_output_functions(fclaw2d_output_functions_t *output_functions)
{
    fclaw2d_output_functions_t *of = output_functions;

    of->f_patch_write_header = &matlab_write_header;
    of->f_patch_write_output = &matlab_write_output;
}

void copy_output_functions(fclaw2d_output_functions_t *old_output_functions,
                           fclaw2d_output_functions_t *new_output_functions)
{
    fclaw2d_output_functions_t* oldof = old_output_functions;
    fclaw2d_output_functions_t* newof = new_output_functions;

    newof->f_patch_write_header = oldof->f_patch_write_header;
    newof->f_patch_write_output = oldof->f_patch_write_output;
}

/* --------------------------------------------------------------
   Functions that the user can replicate
   -------------------------------------------------------------- */

void matlab_write_header(fclaw2d_domain_t* domain, int iframe,int ngrids)
{
    const amr_options_t *gparms = get_domain_parms(domain);
    double time = get_domain_time(domain);

    if (domain->mpirank == 0)
    {
        printf("Matlab output Frame %d  at time %16.8e\n\n",iframe,time);
    }

    // Write out header file containing global information for 'iframe'
    int meqn = gparms->meqn;
    int maux = 0;
    write_tfile_(iframe,time,meqn,ngrids,maux);

    // This opens file 'fort.qXXXX' for replace (where XXXX = <zero padding><iframe>, e.g. 0001,
    // 0010, 0114), and closes the file.
    new_qfile_(iframe);
}


void matlab_write_output(fclaw2d_domain_t *domain, fclaw2d_patch_t *this_patch,
                         int this_block_idx, int this_patch_idx,
                         int iframe,int num,int level)
{
    // In case this is needed by the setaux routine
    set_block_(&this_block_idx);

    /* ----------------------------------------------------------- */
    // Global parameters
    const amr_options_t *gparms = get_domain_parms(domain);
    int mx = gparms->mx;
    int my = gparms->my;
    int mbc = gparms->mbc;
    int meqn = gparms->meqn;

    /* ----------------------------------------------------------- */
    // Patch specific parameters
    ClawPatch *cp = get_clawpatch(this_patch);
    double xlower = cp->xlower();
    double ylower = cp->ylower();
    double dx = cp->dx();
    double dy = cp->dy();

    /* ------------------------------------------------------------ */
    // Pointers needed to pass to Fortran
    double* q = cp->q();

    // Other input arguments
    int maxmx = mx;
    int maxmy = my;

    /* ------------------------------------------------------------- */
    // This opens a file for append.  Now, the style is in the 'clawout' style.
    int matlab_level = level + 1;
    write_qfile_(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,
                 iframe,num,matlab_level,this_block_idx);
}
