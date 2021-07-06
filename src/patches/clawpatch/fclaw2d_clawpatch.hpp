/*
Copyright (c) 2012-2020 Carsten Burstedde, Donna Calhoun
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

#ifndef FCLAW2D_CLAWPATCH_HPP
#define FCLAW2D_CLAWPATCH_HPP

#include <fclaw2d_farraybox.hpp>  /* Needed for FArray boxes */

struct fclaw2d_patch;
struct fclaw2d_global;
class  fclaw2d_metric_patch_t;

class fclaw2d_clawpatch_t
{
public :
    Box dataBox();  /* Box containing data pointer q */
    Box areaBox();  /* Box containing area */
    Box edgeBox();  /* Box containing edge based values */
    Box nodeBox();  /* Box containing nodes */

    /* Solution data */
    int meqn;                   
    FArrayBox griddata;
    FArrayBox griddata_last;
    FArrayBox griddata_save;
    FArrayBox griddata_time_interpolated;
    FArrayBox griderror;

    /* For diagnostics */
    FArrayBox exactsolution;

    int mfields;  /* Number of fields in the rhs */
    FArrayBox rhs;  /* For elliptic problems */

    FArrayBox elliptic_error;  /* For elliptic problems */
    FArrayBox elliptic_soln;  /* For elliptic problems */

    /* Registers for accumulating mismatches at coarse/fine interfaces */
    struct fclaw2d_clawpatch_registers *registers;

    /* Grid info */
    int mx;           
    int my;           
    int mbc;          
    int maux;

    double dx;
    double dy;
    double xlower;
    double ylower;
    double xupper;
    double yupper;

    /* Auxilliary array (used by Clawpack 4.6 and 5.0) */
    FArrayBox aux;
    FArrayBox aux_save;

    /* Mapping and metric info */
    int manifold;    
    int blockno;

    fclaw2d_metric_patch_t *mp;

    /* Extra storage needed by the solver(s) */
    void* solver_data;

    /* User data*/ 
    void* user_data;
};

fclaw2d_clawpatch_t* 
fclaw2d_clawpatch_get_clawpatch(struct fclaw2d_patch* this_patch);

fclaw2d_metric_patch_t* 
fclaw2d_clawpatch_get_metric_patch(struct fclaw2d_patch* this_patch);




#endif /* !FCLAW2D_CLAWPATCH_HPP */
