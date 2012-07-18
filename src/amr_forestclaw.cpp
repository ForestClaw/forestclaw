#include "amr_forestclaw.H"
#include "amr_utils.H"

void amrsetup(fclaw2d_domain_t *domain)
{
    fclaw2d_domain_data_t *dd = P4EST_ALLOC (fclaw2d_domain_data_t, 1);
    domain->user = (void *) dd;

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
            cp->set_mx(domain->mx_leaf);
            cp->set_my(domain->my_leaf);
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
