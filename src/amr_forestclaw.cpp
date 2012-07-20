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

void set_domain_data(fclaw2d_domain_t *domain, global_parms *parms)
{
    fclaw2d_domain_data_t *ddata = P4EST_ALLOC (fclaw2d_domain_data_t, 1);
    domain->user = (void *) ddata;

    ddata->parms = parms;
}

global_parms* get_domain_data(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *ddata;
    ddata = (fclaw2d_domain_data_t *) domain->user;

    return ddata->parms;
}

void set_patch_data(fclaw2d_patch_t *patch, ClawPatch* cp)
{
    fclaw2d_patch_data_t *pdata = P4EST_ALLOC (fclaw2d_patch_data_t, 1);
    patch->user = (void *) pdata;

    pdata->cp = cp;
}

ClawPatch* get_patch_data(fclaw2d_patch_t *patch)
{
    fclaw2d_patch_data_t *pdata;
    pdata = (fclaw2d_patch_data_t *) patch->user;

    return pdata->cp;
}


void amrsetup(fclaw2d_domain_t *domain)
{
    global_parms *gparms = new global_parms();
    cout << "Global parameters " << endl;
    gparms->get_inputParams();
    gparms->print_inputParams();

    set_domain_data(domain,parms);

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = new ClawPatch();

            // Set stuff from p4est

            cp->define(patch->xlower,
                       patch->ylower,
                       patch->xupper,
                       patch->yupper,
                       gparms);
            set_patch_data(patch,cp);
        }
    }
}

void amrinit(fclaw2d_domain_t *domain)
{
    global_parms *gparms = get_domain_data(domain);
    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = get_patch_data(patch);

            cout << "Initializing patch " << j << " on block " << i << endl;
            cp->initialize();
            cp->setAuxArray(gparms->m_maxlevel,gparms->m_refratio,patch->level);

        }
    }
}


void amrrun(fclaw2d_domain_t *domain)
{
    global_parms *parms = get_domain_data(domain);
    double tfinal = parms->m_tfinal;
    cout << "Final time : " << tfinal << endl;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = &domain->blocks[i];
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = &block->patches[j];
            ClawPatch *cp = get_patch_data(patch);

            cout << endl;
            cout << "Clawpatch info for patch number " << j << endl;
            cout << "cp->get_mx() " << cp->get_mx() << endl;
            cout << "cp->get_my() " << cp->get_my() << endl;
            cout << "cp->get_xlower() " << cp->get_xlower() << endl;
            cout << "cp->get_ylower() " << cp->get_ylower() << endl;
            cout << "cp->get_xupper() " << cp->get_xupper() << endl;
            cout << "cp->get_yupper() " << cp->get_yupper() << endl;
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
