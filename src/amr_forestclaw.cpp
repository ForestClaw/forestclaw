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
#include "amr_utils.H"

void amrsetup(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *dd = P4EST_ALLOC (fclaw2d_domain_data_t, 1);
    domain->user = (void *) dd;

    // Todo: Replace this with proper global parameters
    dd->mx_leaf = 8;
    dd->my_leaf = 8;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_data_t *pd = P4EST_ALLOC (fclaw2d_patch_data_t, 1);
            ClawPatch *cp = new ClawPatch();

            patch->user = (void *) pd;
            pd->cp = cp;

            // Set stuff from p4est
            cp->set_mx(dd->mx_leaf);
            cp->set_my(dd->my_leaf);
            cp->set_xlower(patch->xlower);
            cp->set_ylower(patch->ylower);
            cp->set_xupper(patch->xupper);
            cp->set_yupper(patch->yupper);

            // Get the rest of the numerical parameters we need.
            // Todo: Move global parameters into dd
            global_parms parms;
            cp->get_inputParams(parms);
            cp->print_inputParams();
        }
    }
}

void amrrun(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *dd;
    dd = (fclaw2d_domain_data_t *) domain->user;

    global_parms parms;
    ClawPatch *cp1 = new ClawPatch();
    cp1->get_inputParams(parms);

    double tfinal = parms.m_tfinal;
    cout << "Final time : " << tfinal << endl;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            // fclaw2d_patch_t *patch = block->patches + j;
            // fclaw2d_patch_data_t *pd = (fclaw2d_patch_data_t *) patch->user;
            // ClawPatch *cp = pd->cp;
        }
    }
}

void amrreset(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *dd;
    dd = (fclaw2d_domain_data_t *) domain->user;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            fclaw2d_patch_data_t *pd = (fclaw2d_patch_data_t *) patch->user;
            ClawPatch *cp = pd->cp;

            delete cp;
            P4EST_FREE (pd);
            patch->user = NULL;
        }
    }
    P4EST_FREE (dd);
    domain->user = NULL;
}
