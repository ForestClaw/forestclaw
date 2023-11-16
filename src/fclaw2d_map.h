/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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
#endif

#if 0
/* Fix syntax highlighting */
#endif

typedef struct fclaw_global fclaw_global_t;
struct p4est_connectivity;
struct fclaw_global;

struct fclaw_map_context;

/** This prototype matches the Fortran mapc2m functions used in ClawPatch.
 */
typedef void (*fclaw2d_map_c2m_fortran_t) (const double *xc, const double *yc,
                                           double *xp, double *yp,
                                           double *zp);

typedef struct fclaw_map_context fclaw_map_context_t;
typedef struct fclaw2d_map_data fclaw2d_map_data_t;

/** This function is used to query the map for general properties.
 * \param [in] cont     Matching mapping context.
 * \param [in] query_identifier Integer that identifies the query.
 *                      Be sure to use the symbolic constants above.
 * \return              Result of the query.
 */
typedef int (*fclaw_map_query_t) (fclaw_map_context_t * cont,
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
typedef void (*fclaw_map_2d_c2m_t) (fclaw_map_context_t * cont, int blockno,
                                    double xc, double yc,
                                    double *xp, double *yp, double *zp);



/* Covariant and contravariant basis vectors needed for exact solution */
typedef void (*fclaw_map_2d_c2m_basis_t)(fclaw_map_context_t * cont, 
                                         double xc, double yc, 
                                         double *t, double *tinv, 
                                         double *tderivs, int flag);

/* For extruded mesh mappings */
typedef void (*fclaw_map_3d_c2m_t) (fclaw_map_context_t * cont, int blockno,
                                    double xc, double yc,double zc,
                                    double *xp, double *yp, double *zp);


/* Covariant and contravariant basis vectors needed for exact solution */
typedef void (*fclaw_map_3d_c2m_basis_t)(fclaw_map_context_t * cont, 
                                        double xc, double yc, double zc,
                                        double *t, double *tinv, 
                                        double *tderivs, int flag);

/** Destructor for a fclaw2d_map_context.
 */
typedef void (*fclaw_map_destroy_t) (fclaw_map_context_t * cont);

/** Mapping context that is interpreted by its query and c2m members.
 * The callbacks are free to define the meaning of the user_* fields.
 */
struct fclaw_map_context
{
    fclaw_map_query_t       query;

    fclaw_map_2d_c2m_t         mapc2m;
    fclaw_map_2d_c2m_basis_t   basis;

    fclaw_map_3d_c2m_t         mapc2m_3dx;   /* Takes a 2d context */
    fclaw_map_3d_c2m_basis_t   basis_3dx;
    int is_extruded;

    fclaw_map_destroy_t destroy;

    /* Used strictly for 2d mapping */
    int user_int[16];
    double user_double[16];

    /* Create separate data for 2d and 3dx */
    int user_int_3dx[16];
    double user_double_3dx[16];

    double scale[3];
    double shift[3];
    double rotate[9];

    fclaw_map_context_t *brick;
    void *user_data;
};

/**
 * Store the map in the glob
 *
 * @param glob the global context
 * @param map the map to store
 */
void fclaw_map_store (fclaw_global_t* glob,
                      struct fclaw_map_context * map);

/**
 * Get map from the glob
 *
 * @param the global context
 * @return the map
 */
fclaw_map_context_t* fclaw_map_get(fclaw_global_t* glob);



void set_scale(fclaw_map_context_t* cont, const double scale[]);
void set_shift(fclaw_map_context_t* cont, const double shift[]);
void set_rotate(fclaw_map_context_t* cont, const double rotate[]);
void set_default_transform(double scale[],double shift[],double rotate[]);


void scale_map(fclaw_map_context_t* cont,
               double *xp, double *yp, double *zp);
void shift_map(fclaw_map_context_t* cont,
               double *xp, double *yp, double *zp);
void rotate_map(fclaw_map_context_t* cont,
                double *xp, double *yp, double *zp);

#define SET_ROTATION_MATRIX FCLAW_F77_FUNC (set_rotation_matrix,SET_ROTATION_MATRIX)
void SET_ROTATION_MATRIX (const double rot_angles[],double rrot[]);



#define FCLAW_MAP_QUERY FCLAW_F77_FUNC_(fclaw_map_query,FCLAW_MAP_QUERY)


/** Query function for the mapping that can be called from Fortran.
 * \param [in] cont     Mapping context with matching callback functions.
 * \param [in] query_identifier Is passed to the map_query_t function.
 * \param [out] iresult         On return contains result of query.
 */
void FCLAW_MAP_QUERY (fclaw_map_context_t ** cont,
                      const int *query_identifier, int *iresult);



/** Mapping function that can be called from Fortran.
 * \param [in] cont     Mapping context with matching callback functions.
 * \param [in] blockno  Number of the block to be transformed.
 * \param [in] cx       X-coordinate in [block->xlower, block->xupper].
 * \param [in] cy       Y-coordinate in [block->ylower, block->yupper].
 * \param [out] mx      Transformed x-coordinate.
 * \param [out] my      Transformed y-coordinate.
 * \param [out] mz      Transformed z-coordinate.
 */

#define FCLAW_MAP_2D_C2M FCLAW_F77_FUNC_(fclaw_map_2d_c2m,FCLAW_MAP_2D_C2M)
void FCLAW_MAP_2D_C2M (fclaw_map_context_t ** cont, int *blockno,
                       const double *xc, const double *yc,
                       double *xp, double *yp, double *zp);


#define FCLAW_MAP_2D_C2M_BASIS FCLAW_F77_FUNC_(fclaw_map_2d_c2m_basis, \
                                               FCLAW_MAP_2D_C2M_BASIS)

void FCLAW_MAP_2D_C2M_BASIS (fclaw_map_context_t ** cont, 
                             const double *xc, const double *yc,
                             double *t, double *tinv, double *tderivs,
                             int * flag);


#define FCLAW_MAP_3D_C2M FCLAW_F77_FUNC_(fclaw_map_3d_c2m,FCLAW_MAP_3D_C2M)
void FCLAW_MAP_3D_C2M (fclaw_map_context_t ** cont, int *blockno,
                       const double *xc, const double *yc, const double *zc,
                       double *xp, double *yp, double *zp);


#define FCLAW_MAP_3D_C2M_BASIS FCLAW_F77_FUNC_(fclaw_map_3d_c2m_basis, \
                                               FCLAW_MAP_3D_C2M_BASIS)

void FCLAW_MAP_3D_C2M_BASIS (fclaw_map_context_t ** cont, 
                             const double *xc, const double *yc, const double *zc,
                             double *t, double *tinv, double *tderivs,
                             int * flag);


/** Map brick to computational coordinates in [0,1]x[0,1]
 * \param [in] cont     Mapping context with matching callback functions.
 * \param [in] blockno  Number of the block to be transformed.
 * \param [in] cx       X-coordinate in [block->xlower, block->xupper].
 * \param [in] cy       Y-coordinate in [block->ylower, block->yupper].
 * \param [out] mx      Transformed x-coordinate.
 * \param [out] my      Transformed y-coordinate.
 * \param [out] mz      Transformed z-coordinate.
 */

#define FCLAW_MAP_2D_BRICK2C FCLAW_F77_FUNC_(fclaw_map_2d_brick2c, \
                                             FCLAW_MAP_2D_BRICK2C)
void FCLAW_MAP_2D_BRICK2C (fclaw_map_context_t ** cont, int *blockno,
                           const double *xc, const double *yc,
                           double *xp, double *yp, double *zp);


/** Deallocate a mapping context.
 * If the \a destroy member is not NULL, it is called on the context.
 * Otherwise, this function calls FCLAW_FREE (cont).
 * \param [in] cont     Mapping context where the \a destroy member is either
 *                      NULL or a valid function that is then called.
 */
void fclaw_map_destroy (fclaw_map_context_t * cont);

fclaw_map_context_t* fclaw_map_new_nomap_brick(fclaw_map_context_t* brick);
void fclaw_map_2d_c2m_nomap_brick(fclaw_map_context_t * cont, int blockno,
                                  double xc, double yc,
                                  double *xp, double* yp, double *zp);


/* ----------------------------------------------------------------------------------
   Pillowsphere utility
   ---------------------------------------------------------------------------------- */

/* This is called to determine is a map is a pillow sphere or not */
int fclaw_map_pillowsphere(struct fclaw_global* glob);


/* ----------------------------------------------------------------------------------
   Some mapping utility functions
   ---------------------------------------------------------------------------------- */


#define SET_BLOCK FCLAW_F77_FUNC_(set_block,SET_BLOCK)
void SET_BLOCK(const int * a_blockno);


#define FCLAW_MAP_SET_CONTEXT FCLAW_F77_FUNC (fclaw_map_set_context, \
                                              FCLAW_MAP_SET_CONTEXT)
void FCLAW_MAP_SET_CONTEXT (fclaw_map_context_t** a_context);


/* ----------------------------------------------------------------------------------
                                   Headers for mappings
   ---------------------------------------------------------------------------------- */

/* -------------------------------------- No map -------------------------------------- */

fclaw_map_context_t* fclaw_map_new_nomap();


/* --------------------------------- Square mappings ---------------------------------- */

fclaw_map_context_t* fclaw2d_map_new_identity(fclaw_map_context_t *brick);

fclaw_map_context_t* fclaw2d_map_new_cart(fclaw_map_context_t* brick,
                                            const double scale[],
                                            const double shift[]);
  
fclaw_map_context_t* fclaw2d_map_new_fivepatch(const double scale[],
                                                 const double shift[],
                                                 const double alpha);

fclaw_map_context_t* fclaw2d_map_new_squareddisk(const double scale[],
                                                   const double shift[],
                                                   const double alpha);
  
fclaw_map_context_t* fclaw2d_map_new_bilinear(fclaw_map_context_t *brick,
                                                const double scale[],
                                                const double shift[],
                                                const double center[]);

/* ---------------------------------- Disk mappings ----------------------------------- */

fclaw_map_context_t* fclaw2d_map_new_pillowdisk(const double scale[],
                                                  const double shift[],
                                                  const double rotate[]);

fclaw_map_context_t* fclaw2d_map_new_pillowdisk5(const double scale[],
                                                   const double shift[],
                                                   const double rotate[],
                                                   const double alpha);

/* --------------------------------- Annulus mapping ---------------------------------- */
fclaw_map_context_t *
    fclaw2d_map_new_annulus (fclaw_map_context_t* brick,
                             const double scale[],
                             const double rotate[],
                             const double alpha,
                             const double theta[]);


/* --------------------------------- Latlong mapping ---------------------------------- */
fclaw_map_context_t *
    fclaw2d_map_new_latlong (fclaw_map_context_t* brick,
                             const double scale[],
                             const double rotate[],
                             const double lat[],
                             const double longitude[],
                             const int a, const int b);

/* --------------------------------- Hemisphere mappings ------------------------------ */

fclaw_map_context_t* fclaw2d_map_new_pillowsphere5(const double scale[],
                                                     const double rotate[],
                                                     const double alpha);

/* --------------------------------- Sphere mappings ---------------------------------- */

fclaw_map_context_t* fclaw2d_map_new_pillowsphere(const double scale[],
                                                    const double rotate[]);

fclaw_map_context_t * fclaw2d_map_new_cubedsphere (const double scale[],
                                                     const double rotate[]);

/* --------------------------------- Torus mappings ---------------------------------- */
fclaw_map_context_t *
    fclaw2d_map_new_torus (fclaw_map_context_t* brick,
                           const double scale[],
                           const double rotate[],
                           const double alpha,
                           const double beta);





/* ----------------------------------------------------------------------------------
   Some generic fortran mappings.  Users can call these by setting up a
   'fclaw2d_map_<name>_new() function.
   ---------------------------------------------------------------------------------- */

/* Single block mappings */
#define FCLAW_MAP_2D_C2M_IDENTITY FCLAW_F77_FUNC (fclaw_map_2d_c2m_identity,FCLAW_MAP_2D_C2M_IDENTITY)
void FCLAW_MAP_2D_C2M_IDENTITY (int* blockno, double *xc, double *yc,
                      double *xp, double *yp, double *zp);

/* Single block mappings */
#define FCLAW_MAP_2D_C2M_CART FCLAW_F77_FUNC (fclaw_map_2d_c2m_cart,FCLAW_MAP_2D_C2M_CART)
void FCLAW_MAP_2D_C2M_CART (int* blockno, double *xc, double *yc,
                  double *xp, double *yp, double *zp);


#define FCLAW_MAP_2D_C2M_PILLOWDISK FCLAW_F77_FUNC (fclaw_map_2d_c2m_pillowdisk,FCLAW_MAP_2D_C2M_PILLOWDISK)
void FCLAW_MAP_2D_C2M_PILLOWDISK (int* blockno, double *xc, double *yc,
                        double *xp, double *yp, double *zp);

#define FCLAW_MAP_2D_C2M_PILLOWDISK5 FCLAW_F77_FUNC (fclaw_map_2d_c2m_pillowdisk5,FCLAW_MAP_2D_C2M_PILLOWDISK5)
void FCLAW_MAP_2D_C2M_PILLOWDISK5 (int* blockno, double *xc, double *yc,
                         double *xp, double *yp, double *zp, double *alpha);

/* multi-block mappings */
#define FCLAW_MAP_2D_C2M_SQUAREDDISK FCLAW_F77_FUNC (fclaw_map_2d_c2m_squareddisk,FCLAW_MAP_2D_C2M_SQUAREDDISK)
void FCLAW_MAP_2D_C2M_SQUAREDDISK (int *blockno, double *xc, double *yc,
                         double *xp, double *yp, double *zp, double *alpha);

#define FCLAW_MAP_2D_C2M_FIVEPATCH FCLAW_F77_FUNC (fclaw_map_2d_c2m_fivepatch,FCLAW_MAP_2D_C2M_FIVEPATCH)
void FCLAW_MAP_2D_C2M_FIVEPATCH (int* blockno, double *xc, double *yc,
                       double *xp, double *yp, double *zp,double *alpha);

#define FCLAW_MAP_2D_C2M_CUBEDSPHERE FCLAW_F77_FUNC (fclaw_map_2d_c2m_cubedsphere,FCLAW_MAP_2D_C2M_CUBEDSPHERE)
void FCLAW_MAP_2D_C2M_CUBEDSPHERE (int* blockno, double *xc, double *yc,
                         double *xp, double *yp, double *zp);

#define FCLAW_MAP_2D_C2M_PILLOWSPHERE FCLAW_F77_FUNC (fclaw_map_2d_c2m_pillowsphere,FCLAW_MAP_2D_C2M_PILLOWSPHERE)
void FCLAW_MAP_2D_C2M_PILLOWSPHERE (int* blockno, double *xc, double *yc,
                          double *xp, double *yp, double *zp);

#define FCLAW_MAP_2D_C2M_TORUS FCLAW_F77_FUNC (fclaw_map_2d_c2m_torus,FCLAW_MAP_2D_C2M_TORUS)
void FCLAW_MAP_2D_C2M_TORUS (double *xc, double *yc, double *xp, double *yp, double *zp, 
                   double* alpha, double* beta);

#define FCLAW_MAP_2D_C2M_TWISTED_TORUS FCLAW_F77_FUNC (fclaw_map_2d_c2m_twisted_torus,FCLAW_MAP_2D_C2M_TWISTED_TORUS)
void FCLAW_MAP_2D_C2M_TWISTED_TORUS (int* blockno, double *xc, double *yc,
                   double *xp, double *yp, double *zp, double* alpha);

#define FCLAW_MAP_2D_C2M_BRICK FCLAW_F77_FUNC (fclaw_map_2d_c2m_brick,FCLAW_MAP_2D_C2M_BRICK)
void FCLAW_MAP_2D_C2M_BRICK (int* blockno, double *xc, double *yc,
                   double *xp, double *yp, double *zp, int *mi, int *mj);

#define FCLAW_MAP_2D_C2M_LATLONG FCLAW_F77_FUNC (fclaw_map_2d_c2m_latlong,FCLAW_MAP_2D_C2M_LATLONG)
void FCLAW_MAP_2D_C2M_LATLONG (int* blockno, double *xc, double *yc,
                   double *xp, double *yp, double *zp);

#define FCLAW_MAP_2D_C2M_ANNULUS FCLAW_F77_FUNC (fclaw_map_2d_c2m_annulus,FCLAW_MAP_2D_C2M_ANNULUS)
void FCLAW_MAP_2D_C2M_ANNULUS (int* blockno, double *xc, double *yc,
                     double *xp, double *yp, double *zp, double *alpha,
                     double *theta);

/* ---------------------------------------------------------------------------------- */

#ifdef __cplusplus
}
#endif

#endif
