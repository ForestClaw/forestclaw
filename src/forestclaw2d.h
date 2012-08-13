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
typedef struct fclaw2d_block fclaw2d_block_t;
typedef struct fclaw2d_patch fclaw2d_patch_t;

typedef void (*fclaw2d_mapc2m_t) (const double xyc[2], double xyzp[3],
				  fclaw2d_domain_t *domain, void *user);

/* The domain structure gives a processor local view of the grid hierarchy */

struct fclaw2d_patch
{
  int			level;		/* 0 is root, increases if refined */
  double		xlower, xupper;
  double		ylower, yupper;
  fclaw2d_patch_t	*next;		/* next patch same level same block */
  void			*user;
};

struct fclaw2d_block
{
  double		xlower, xupper;
  double		ylower, yupper;
  fclaw2d_mapc2m_t	mapc2m;
  void			*mapc2m_user;
  int			mthbc[4];	/* >0 for physical bc types */
  int			num_patches;	/* proc-local patches in this block */
  int			num_patches_before;	/* in all previous blocks */
  fclaw2d_patch_t	*patches;		/* allocated storage */
  fclaw2d_patch_t	*patchbylevel[P4EST_MAXLEVEL + 1];	/* pointers */
  void			*user;
};

struct fclaw2d_domain
{
  int			mpisize, mpirank;	/* MPI variables */
  int			num_patches_all;	/* sum over all blocks */
  int			maxlevel_all;		/* maximum over all blocks */
  int			num_blocks;
  fclaw2d_block_t	*blocks;		/* allocated storage */
  int			*patch_to_block;	/* allocated storage */
  p4est_wrap_t		*pp;
  void			*user;
};

typedef void (*fclaw2d_patch_callback_t)
			(fclaw2d_domain_t *domain, fclaw2d_patch_t *patch,
			 int blockno, int patchno, void *user);

void fclaw2d_domain_iterate_level (fclaw2d_domain_t *domain, int level,
                                   fclaw2d_patch_callback_t pcb, void *user);

void fclaw2d_domain_iterate_patches (fclaw2d_domain_t *domain,
                                     fclaw2d_patch_callback_t pcb, void *user);



/** Determine boundary type >0 from fclaw2d_block_t, or 0 for neighbor patches.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in,out] boundaries	Domain boundaries as present in block->mthbc.
 *			The order is left, right, bottom, top.
 * \return		True if at least one patch face is on a boundary.
 */
int			fclaw2d_patch_boundary_type (fclaw2d_domain_t *domain,
						int blockno, int patchno,
						int boundaries[4]);

typedef enum fclaw2d_face_neighbor
{
  FCLAW2D_FACE_NEIGHBOR_BOUNDARY,
  FCLAW2D_FACE_NEIGHBOR_HALFSIZE,
  FCLAW2D_FACE_NEIGHBOR_SAMESIZE,
  FCLAW2D_FACE_NEIGHBOR_DOUBLESIZE
}
fclaw2d_face_neighbor_t;

/** Determine neighbor patch(es) and orientation across a given face.
 * \param [in] domain	Valid domain structure.
 * \param [in] blockno	Number of the block within the domain.
 * \param [in] patchno	Number of the patch within the block.
 * \param [in] faceno	Number of the patch face: left, right, bottom, top.
 * \param [in,out] rproc	Processor number of neighbor patches.
 * \param [in,out] rblockno	Neighbor block number for up to 2 neighbors.
 * \param [in,out] rpatchno	Neighbor patch number for up to 2 neighbors.
 * \param [in,out] rfaceno	Neighbor face number and orientation.
 * \return			The Type of face neighbor connection.
 */
fclaw2d_face_neighbor_t	fclaw2d_patch_face_neighbors (fclaw2d_domain_t *domain,
                                                      int blockno, int patchno, int faceno,
                                                      int rproc[2], int rblockno[2],
                                                      int rpatchno[2], int *rfaceno);

/** Mark a patch for refinement.
 */
void			fclaw2d_patch_mark_refine (fclaw2d_domain_t *domain,
				int blockno, int patchno);

/** Mark a patch for coarsening.
 * Coarsening will only happen if all sibling patches are marked as well.
 */
void			fclaw2d_patch_mark_coarsen (fclaw2d_domain_t *domain,
				int blockno, int patchno);

#ifdef __cplusplus
}
#endif

#endif
