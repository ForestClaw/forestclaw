#include "ClawPatch.H"
#include "forestclaw2d.h"
#include "amr_utils.H"

void amrsetup(fclaw2d_domain_t *domain)
{

    global_parms parms;

    inputparms_(parms.m_initial_dt,
                parms.m_tfinal,
                parms.m_max_cfl,
                parms.m_nout,
                parms.m_src_term,
                parms.m_mcapa,
                parms.m_maux,
                parms.m_meqn,
                parms.m_mwaves,
                parms.m_maxmwaves,
                parms.m_mthlim,
                parms.m_mbc,
                parms.m_mthbc,
                parms.m_order);

    for(int i = 0; i < domain->num_blocks; i++)
    {
        fclaw2d_block_t *block = domain->blocks + i;
        for(int j = 0; j < block->num_patches; j++)
        {
            fclaw2d_patch_t *patch = block->patches + j;
            ClawPatch *cp = new ClawPatch();
            patch->user = cp;

            cp->set_mx(domain->mx_leaf);
            cp->set_my(domain->my_leaf);

            // cp->set_order(parms.m_order);
            cp->set_src(parms.m_src_term);
            cp->set_mcapa(parms.m_mcapa);
            cp->set_maux(parms.m_maux);
            cp->set_mbc(parms.m_mbc);
            cp->set_meqn(parms.m_meqn);
            cp->set_mwaves(parms.m_mwaves);
            // cp->set_mthlim(parms.m_mthlim);
            // cp->set_mthbc(parms.m_mthbc);
            cp->set_initial_dt(parms.m_initial_dt);
            cp->set_xlower(patch->xlower);
            cp->set_ylower(patch->ylower);

            cp->set_xupper(patch->xupper);
            cp->set_yupper(patch->yupper);

            cp->set_max_cfl(parms.m_max_cfl);
        }
    }
}
