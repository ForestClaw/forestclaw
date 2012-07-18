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

#ifndef FORESTCLAW2D_H
#define FORESTCLAW2D_H

#include <p4est_wrap.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const double fclaw2d_root_len;
extern const double fclaw2d_smallest_h;

typedef struct fclaw2d_domain fclaw2d_domain_t;

typedef void (*fclaw2d_mapc2m_t) (const double xyc[2], double xyzp[P4EST_DIM],
				  fclaw2d_domain_t *domain, void *user);

typedef struct fclaw2d_patch
{
  int			level;		/* 0 is root, increases if refined */
  double		xlower, xupper;
  double		ylower, yupper;
  void			*user;
}
fclaw2d_patch_t;

typedef struct fclaw2d_block
{
  double		xlower, xupper;
  double		ylower, yupper;
  fclaw2d_mapc2m_t	mapc2m;
  void			*mapc2m_user;
  int                   mthbc[P4EST_FACES];	/* >0 for physical bc types */
  int			num_patches;
  fclaw2d_patch_t	*patches;
  void			*user;
}
fclaw2d_block_t;

struct fclaw2d_domain
{
  int			num_blocks;
  fclaw2d_block_t	*blocks;
  p4est_wrap_t          *pp;
  void			*user;
};

/** Determine boundary type >0 from fclaw2d_block_t, or 0 for neighbor patches.
 * \param [in]		number of the block within the domain.
 * \param [in]		number of the patch within the block.
 * \param [in,out]	boundary types as present in fclaw2d_block_t.
 *			The order is left, right, bottom, top.
 * \return		True if at least one patch face is on a boundary.
 */
int			fclaw2d_patch_boundary_type (fclaw2d_domain_t *domain,
						int blockno, int patchno,
						int boundaries[P4EST_FACES]);

#ifdef __cplusplus
}
#endif

#endif
