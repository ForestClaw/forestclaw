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

#ifndef FCLAW2D_MAP_H
#define FCLAW2D_MAP_H

#include <fclaw_base.h>
#include <fclaw2d_map_query_defs.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif




/** This prototype matches the Fortran mapc2m functions used in ClawPatch.
 */
typedef void (*fclaw2d_map_c2m_fortran_t) (const double *xc, const double *yc,
                                           double *xp, double *yp,
                                           double *zp);

#define FCLAW2D_MAP_QUERY FCLAW_F77_FUNC_(fclaw2d_map_query,FCLAW2D_MAP_QUERY)


typedef struct fclaw2d_map_context fclaw2d_map_context_t;

/** This function is used to query the map for general properties.
 * \param [in] cont     Matching mapping context.
 * \param [in] query_identifier Integer that identifies the query.
 *                      Be sure to use the symbolic constants above.
 * \return              Result of the query.
 */
typedef int (*fclaw2d_map_query_t) (fclaw2d_map_context_t * cont,
                                    int query_identifier);

/** This function performs the coordinate transformation.
 * \param [in] cont     Matching mapping context.
 * \param [in] blockno  Number of the block to be transformed.
 * \param [in] xc       X-coordinate in [block->xlower, block->xupper].
 * \param [in] yc       Y-coordinate in [block->ylower, block->yupper].
 * \param [out] xp      Transformed x-coordinate.
 * \param [out] yp      Transformed y-coordinate.
 * \param [out] zp      Transformed z-coordinate.
 */
typedef void (*fclaw2d_map_c2m_t) (fclaw2d_map_context_t * cont, int blockno,
                                   double xc, double yc,
                                   double *xp, double *yp, double *zp);

/** Destructor for a fclaw2d_map_context.
 */
typedef void (*fclaw2d_map_destroy_t) (fclaw2d_map_context_t * cont);

/** Mapping context that is interpreted by its query and c2m members.
 * The callbacks are free to define the meaning of the user_* fields.
 */
struct fclaw2d_map_context
{
    fclaw2d_map_query_t query;
    fclaw2d_map_c2m_t mapc2m;
    fclaw2d_map_destroy_t destroy;
    int user_int[16];
    double user_double[16];
    void *user_data;
};


/** Query function for the mapping that can be called from Fortran.
 * \param [in] cont     Mapping context with matching callback functions.
 * \param [in] query_identifier Is passed to the map_query_t function.
 * \param [out] iresult         On return contains result of query.
 */


void FCLAW2D_MAP_QUERY (fclaw2d_map_context_t ** cont,
                        const int *query_identifier, int *iresult);


#define FCLAW2D_MAP_C2M FCLAW_F77_FUNC_(fclaw2d_map_c2m,FCLAW2D_MAP_C2M)

/** Mapping function that can be called from Fortran.
 * \param [in] cont     Mapping context with matching callback functions.
 * \param [in] blockno  Number of the block to be transformed.
 * \param [in] cx       X-coordinate in [block->xlower, block->xupper].
 * \param [in] cy       Y-coordinate in [block->ylower, block->yupper].
 * \param [out] mx      Transformed x-coordinate.
 * \param [out] my      Transformed y-coordinate.
 * \param [out] mz      Transformed z-coordinate.
 */
void FCLAW2D_MAP_C2M (fclaw2d_map_context_t ** cont, int *blockno,
                      const double *xc, const double *yc,
                      double *xp, double *yp, double *zp);

/** Deallocate a mapping context.
 * If the \a destroy member is not NULL, it is called on the context.
 * Otherwise, this function calls FCLAW_FREE (cont).
 * \param [in] cont     Mapping context where the \a destroy member is either
 *                      NULL or a valid function that is then called.
 */
void fclaw2d_map_destroy (fclaw2d_map_context_t * cont);


/* ----------------------------------------------------------------------------------
   New maps (torus, cubedsphere, disk) and a utility function for calling maps
   defined in fortran.
   ---------------------------------------------------------------------------------- */

/** Create a torus mapping for one block with [0, 1]^2 (for now).
 * \param [in] R1       Large radius of the torus.
 * \param [in] R2       Small radius of the torus.
 * \return              Mapping context.
 */
fclaw2d_map_context_t *fclaw2d_map_new_torus (double R1, double R2);

/** Create a cubed sphere mapping from six trees.
 * \param [in] R        Radius of the cubed sphere surface.
 * \return              Mapping context.
 */
fclaw2d_map_context_t *fclaw2d_map_new_csphere (double R);

/** Create a planar spherical disk mapping from five trees.
 * It is composed of a center square and one deformed patch on either side.
 * \param [in] R1       Outer radius of the disk.
 * \param [in] R2       Radius at any corner of the inside square.
 * \return              Mapping context.
 */
fclaw2d_map_context_t *fclaw2d_map_new_disk (double R1, double R2);

/** Create a rotated, scaled pillow disk;
 * It is composed of a center square and one deformed patch on either side.
 * \param [in] rotate   (theta,phi) rotation angles
 * \param [in] scale    Scale the unit cube (e.g. set the radius)
 * \return              Mapping context.
 */

fclaw2d_map_context_t * fclaw2d_map_new_pillowdisk (double rotate[], double scale);

/** Create a rotated, scaled squared disk;
 * It is composed of a center square and one deformed patch on either side.
 * \param [in] rotate   (theta,phi) rotation angles
 * \param [in] scale    Scale the unit cube (e.g. set the radius)
 * \return              Mapping context.
 */
fclaw2d_map_context_t * fclaw2d_map_new_squareddisk (double rotate[], double scale);



/** Create a rotated, scaled cubed sphere. Uses same cubed sphere map as above.
 * It is composed of a center square and one deformed patch on either side.
 * \param [in] rotate   (theta,phi) rotation angles
 * \param [in] scale    Scale the unit cube (e.g. set the radius)
 * \return              Mapping context.
 */
fclaw2d_map_context_t * fclaw2d_map_new_cubedsphere (double rotate[], double scale);


/** Create a rotated, scaled pillow sphere.
 * It is composed of a center square and one deformed patch on either side.
 * \param [in] rotate   (theta,phi) rotation angles
 * \param [in] scale    Scale the unit cube (e.g. set the radius)
 * \return              Mapping context.
 */
fclaw2d_map_context_t * fclaw2d_map_new_pillowsphere (double rotate[], double scale);


/** Create a mapping context for any number of blocks using a Fortran mapc2m.
 * \param [in] mapc2m   Address of the Fortran mapping function.
 *                      It expects the block number in a Clawpatch COMMON.
 * \param [in] query_results    Results for the queries defined above.
 *                      Be sure to assign them by the symbolic constants
 *                      defined above since these may change in the future.
 * \return              Mapping context.
 */
fclaw2d_map_context_t *fclaw2d_map_new_fortran (fclaw2d_map_c2m_fortran_t
                                                mapc2m,
                                                const int
                                                query_results
                                                [FCLAW2D_MAP_QUERY_LAST]);

/* ----------------------------------------------------------------------------------
   Mapping routines (from clawpack_fort.H)
   ---------------------------------------------------------------------------------- */

#define ISPILLOWSPHERE FCLAW_F77_FUNC_(ispillowsphere,ISPILLOWSPHERE)
int ISPILLOWSPHERE();

#define IS_CUBEDSPHERE FCLAW_F77_FUNC_(iscubedsphere,ISCUBEDSPHERE)
int ISCUBEDSPHERE();

#define ISPILLOWDISK FCLAW_F77_FUNC_(ispillowdisk,ISPILLOWDISK)
int ISPILLOWDISK();

#define ISFLAT FCLAW_F77_FUNC_(isflat,ISFLAT)
int ISFLAT();

#define ISSPHERE FCLAW_F77_FUNC_(issphere,ISSPHERE)
int ISSPHERE();

#define SET_BLOCK FCLAW_F77_FUNC_(set_block,SET_BLOCK)
void SET_BLOCK(const int * a_blockno);

#define SET_CONTEXT FCLAW_F77_FUNC (set_context,SET_CONTEXT)
void SET_CONTEXT (fclaw2d_map_context_t** a_context);

#define SETUP_MAPPEDGRID FCLAW_F77_FUNC (setup_mappedgrid,SETUP_MAPPEDGRID)
void SETUP_MAPPEDGRID(double rot_angle[], double* scale);

#define SCALE_MAP FCLAW_F77_FUNC (scale_map,SCALE_MAP)
void SCALE_MAP (double *xp, double *yp, double *zp);

#define ROTATE_MAP FCLAW_F77_FUNC (rotate_map,ROTATE_MAP)
void ROTATE_MAP (double *xp, double *yp, double *zp);

#define MAPC2M_PILLOWDISK FCLAW_F77_FUNC (mapc2m_pillowdisk,MAPC2M_PILLOWDISK)
void MAPC2M_PILLOWDISK (*xc, double *yc, double *xp, double *yp, double *zp);

#define MAPC2M_CUBEDSPHERE FCLAW_F77_FUNC (mapc2m_cubedsphere,MAPC2M_CUBEDSPHERE)
void MAPC2M_CUBEDSPHERE (double *xc, double *yc, double *xp, double *yp, double *zp);

#define MAPC2M_PILLOWSPHERE FCLAW_F77_FUNC (mapc2m_pillowsphere,MAPC2M_PILLOWSPHERE)
void MAPC2M_PILLOWSPHERE (double *xc, double *yc, double *xp, double *yp, double *zp);
/* ---------------------------------------------------------------------------------- */



#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
