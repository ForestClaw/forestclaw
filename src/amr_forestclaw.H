#ifndef AMR_FORESTCLAW_H
#define AMR_FORESTCLAW_H

#include "ClawPatch.H"
#include "forestclaw2d.h"

typedef struct fclaw2d_domain_data
{
  // int		mx_leaf, my_leaf;
}
fclaw2d_domain_data_t;

typedef struct fclaw2d_patch_data
{
  ClawPatch	*cp;
}
fclaw2d_patch_data_t;

void amrsetup(fclaw2d_domain_t *domain);
void amrrun(fclaw2d_domain_t *domain);
void amrreset(fclaw2d_domain_t *domain);

#endif
