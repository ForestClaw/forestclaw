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

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

/** This prototype matches the Fortran mapc2m functions used in ForestClaw.
 */
typedef void (*fclaw2d_map_c2m_fortran_t) (const double *cx, const double *cy,
                                           double *mx, double *my,
                                           double *mz);

/* These integers are meant to be passed in query_identifier of map_query_t.
 * One of these four types is generally chosen.
 * 1 maps into R^2, 2 and 3 map into R^3.
 */
#define FCLAW2D_MAP_QUERY_IS_USED         0     /* is the map used at all? */
#define FCLAW2D_MAP_QUERY_IS_SCALEDSHIFT  1     /* x_i -> a_i x_i + b_i */
#define FCLAW2D_MAP_QUERY_IS_AFFINE       2     /* x -> A x + b */
#define FCLAW2D_MAP_QUERY_IS_NONLINEAR    3     /* x -> F(x) */

/* Query additional properties of mappings.
 * These properties can be used to implement shortcuts in the numerical code.
 */
#define FCLAW2D_MAP_QUERY_IS_GRAPH        4     /* (x,y) -> (x,y,f(x,y)) */
#define FCLAW2D_MAP_QUERY_LAST            5     /* #"official" queries. */

typedef struct fclaw2d_map_context fclaw2d_map_context_t;

/** This function is used to query the map for general properties.
 * \param [in] cont     Matching mapping context.
 * \param [in] query_identifier Integer that identifies the query.
 * \return              Result of the query.
 */
typedef int (*fclaw2d_map_query_t) (fclaw2d_map_context_t * cont,
                                    int query_identifier);

/** This function performs the coordinate transformation.
 * \param [in] cont     Matching mapping context.
 * \param [in] blockno  Number of the block to be transformed.
 * \param [in] cx       X-coordinate in [block->xlower, block->xupper].
 * \param [in] cy       Y-coordinate in [block->ylower, block->yupper].
 * \param [out] mx      Transformed x-coordinate.
 * \param [out] my      Transformed y-coordinate.
 * \param [out] mz      Transformed z-coordinate.
 */
typedef void (*fclaw2d_map_c2m_t) (fclaw2d_map_context_t * cont, int blockno,
                                   double cx, double cy,
                                   double *mx, double *my, double *mz);

/** Mapping context that is interpreted by its query and c2m members.
 * The callbacks are free to define the meaning of the user_* fields.
 */
struct fclaw2d_map_context
{
    unsigned magic;
    fclaw2d_map_query_t query;
    fclaw2d_map_c2m_t mapc2m;
    void *user_data;
    int user_int[16];
    double user_double[16];
};

/** Query function for the mapping that can be called from Fortran.
 * \param [in] cont     Mapping context with matching callback functions.
 * \param [in] query_identifier Is passed to the map_query_t function.
 * \param [out] iresult         On return contains result of query.
 */
void fclaw2d_map_query_ (fclaw2d_map_context_t * cont,
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
void fclaw2d_map_c2m_ (fclaw2d_map_context_t * cont, int *blockno,
                       const double *cx, const double *cy,
                       double *mx, double *my, double *mz);

/** Create a torus mapping for one block with [0, 1]^2 (for now).
 * \param [in] R1       Large radius of the torus.
 * \param [in] R2       Small radius of the torus.
 * \return              Mapping context.
 *                      Must be destroyed by fclaw2d_map_destroy_torus.
 */
fclaw2d_map_context_t *fclaw2d_map_new_torus (double R1, double R2);
void fclaw2d_map_destroy_torus (fclaw2d_map_context_t * cont);

/** Create a mapping context for any number of blocks using a Fortran mapc2m.
 * \param [in] mapc2m   Address of the Fortran mapping function.
 *                      It expects the block number in a Clawpatch COMMON.
 * \param [in] query_results    Results for the queries defined above.
 * \return              Mapping context.
 *                      Must be destroyed by fclaw2d_map_destroy_fortran.
 */
fclaw2d_map_context_t *fclaw2d_map_new_fortran (fclaw2d_map_c2m_fortran_t
                                                mapc2m,
                                                const int
                                                query_results
                                                [FCLAW2D_MAP_QUERY_LAST]);
void fclaw2d_map_destroy_fortran (fclaw2d_map_context_t * cont);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif
