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

#ifndef FCLAW2D_CONVENIENCE_H
#define FCLAW2D_CONVENIENCE_H

#include <forestclaw2d.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}                               /* need this because indent is dumb */
#endif
#endif

fclaw2d_domain_t *fclaw2d_domain_new_unitsquare (sc_MPI_Comm mpicomm,
                                                 int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_torus (sc_MPI_Comm mpicomm,
                                            int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_twosphere (sc_MPI_Comm mpicomm,
                                                int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_cubedsphere (sc_MPI_Comm mpicomm,
                                                  int initial_level);
fclaw2d_domain_t *fclaw2d_domain_new_disk (sc_MPI_Comm mpicomm,
                                           int initial_level);

void fclaw2d_domain_destroy (fclaw2d_domain_t * domain);

/** Create a new domain based on refine and coarsen marks set previously.
 * All refine and coarsen markers are cancelled when this function is done.
 * \param [in,out] domain       Current domain with set adaptation markers.
 *                              It stays alive because it is needed to
 *                              project numerical values to the adapted domain.
 *                              If adapted, no queries are allowed afterwards.
 * \return                      Adapted domain if refinement occurred, or NULL.
 *                              The return status is identical across all ranks.
 */
fclaw2d_domain_t *fclaw2d_domain_adapt (fclaw2d_domain_t * domain);

/** Create a repartitioned domain after fclaw2d_domain_adapt returned non-NULL.
 * All refine and coarsen markers are cancelled when this function is done.
 * \param [in,out] domain       Current domain that was adapted previously.
 *                              It stays alive because it is needed to
 *                              transfer numerical values to the new partition.
 *                              If partitioned, no queries allowed afterwards.
 * \param [in] weight_exponent  The patches are weighted with an integer factor
 *                              2 ** (level * exponent).  If the exponent is 0,
 *                              all patches have equal weight.  If it is 1,
 *                              smaller patches are considered more expensive
 *                              by a factor two per level increase.
 * \return                      Partitioned domain if different, or NULL.
 *                              The return status is identical across all ranks.
 */
fclaw2d_domain_t *fclaw2d_domain_partition (fclaw2d_domain_t * domain,
                                            int weight_exponent);

/** Clean up after fclaw2d_domain_partition returned non-NULL.
 * \param [in,out] domain       Current domain that was partitioned.
 */
void fclaw2d_domain_complete (fclaw2d_domain_t * domain);

/** Print patch number by level on all processors */
void fclaw2d_domain_list_levels (fclaw2d_domain_t * domain, int log_priority);

/** Print face neighbor status for each face */
void fclaw2d_domain_list_neighbors (fclaw2d_domain_t * domain,
                                    int log_priority);

/** Print information on adapted patches */
void fclaw2d_domain_list_adapted (fclaw2d_domain_t * old_domain,
                                  fclaw2d_domain_t * new_domain,
                                  int log_priority);

#ifdef __cplusplus
#if 0
{                               /* need this because indent is dumb */
#endif
}
#endif

#endif /* !FCLAW2D_CONVENIENCE_H */
