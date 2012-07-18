#include "ClawPatch.H"
#include "forestclaw2d.h"
#include "amr_utils.H"

void amrsetup(pfclaw_domain_t *domain)
{
    for(int i = 0; i < domain->num_blocks; i++)
    {
        pfclaw_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            pfclaw_patch_t *patch = block->patches + j;
            ClawPatch *cp = new ClawPatch();
            patch->user = cp;

            cp->set_mx(domain->mx);
            cp->set_my(domain->my);
        }
    }
}
