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

#include <fclaw2d_to_3d.h>
#include "fclaw2d_patch.c"

void fclaw3d_patch_set_edge_type(fclaw2d_patch_t *patch,int iedge,
								   fclaw2d_patch_relation_t edge_type)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	pdata->edge_neighbors[iedge] = edge_type;
}

fclaw2d_patch_relation_t fclaw3d_patch_get_edge_type(fclaw2d_patch_t* patch,
													   int iedge)
{
	fclaw2d_patch_data_t *pdata = get_patch_data(patch);
	FCLAW_ASSERT(pdata->neighbors_set != 0);
	FCLAW_ASSERT(0 <= iedge && iedge < FCLAW3D_NUMEDGES);
	return pdata->edge_neighbors[iedge];
}