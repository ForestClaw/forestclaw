#include "ClawPatch.H"
#include "forestclaw2d.h"
#include "amr_utils.H"
#include "amr_forestclaw.H"

void amrsetup(fclaw2d_domain_t *domain)
{

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            ClawPatch *cp = new ClawPatch();
            patch->user = cp;

            global_parms *parms = new global_parms;

            // Set stuff from p4est
            cp->set_mx(domain->mx_leaf);
            cp->set_my(domain->my_leaf);
            cp->set_xlower(patch->xlower);
            cp->set_ylower(patch->ylower);
            cp->set_xupper(patch->xupper);
            cp->set_yupper(patch->yupper);

            // Get the rest of the numerical parameters we need.
            cp->get_inputParams(*parms);
            cp->print_inputParams();
            block->user = cp;
        }
    }
}

void amrrun(fclaw2d_domain_t *domain)
{
    // Get first block so we can get global parameters
    fclaw2d_block_t *block = domain->blocks;

    global_parms *parms = (global_parms*) block->user;

    double tfinal = parms->m_tfinal;

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            ClawPatch *cp = (ClawPatch*) patch->user;

        }
    }
}
