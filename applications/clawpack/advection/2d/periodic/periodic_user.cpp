/*
Copyright (c) 2012-2025 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include "periodic_user.h"
#include "math.h"

static
void periodic_problem_setup(fclaw_global_t* glob)
{
    const user_options_t* user = periodic_get_options(glob);

    const fclaw_options_t * fclaw_opt = fclaw_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob.data","w");
        fprintf(f,"%-24.4f %s\n",user->uvel,"\% u-velocity");
        fprintf(f,"%-24.4f %s\n",user->vvel,"\% v-velocity");
        fprintf(f,"%-24d %s\n",fclaw_opt->moving_gauges,"\% moving_gauges");
        fclose(f);
    }

    /* Make sure that node 0 has written 'setprob.data' before proceeding */
    fclaw_domain_barrier(glob->domain);

    SETPROB();  /* Reads file created above */
}

void periodic_gauge_move(fclaw_global_t* glob, fclaw_gauge_t *g,
                         double t, double dt)
{
    /* These will reposition all gauges, even if they are currently not 
       in the domain */

    double xc,yc;
    if (g->num == 0)
    {
        /* Prescribed velocity : gauge travels in a circle */
        xc = 0.5*cos(M_PI*t);
        yc = 0.5*sin(M_PI*t);        
    }
    else
    {
        /* Gauge travels in straight line.  May start outside of the domain. */
        int num, dim;
        double zc,t1,t2;            
        fclaw_gauges_get_data(glob,g,&num,&dim,&xc,&yc,&zc,&t1,&t2);

        const user_options_t* user = periodic_get_options(glob);

        /* With this velocity, the gauges may leave the domain before the 
           simulation is done. */
        xc += dt*user->uvel;
        yc += dt*user->vvel;
    }

    fclaw_gauges_set_position(glob, g, xc, yc, g->zc);
}


void periodic_link_solvers(fclaw_global_t *glob)
{
    fclaw_vtable_t *vt = fclaw_vt(glob);
    vt->problem_setup = &periodic_problem_setup;  /* Version-independent */

    const user_options_t* user = periodic_get_options(glob);
    if (user->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t *clawpack46_vt = fc2d_clawpack46_vt(glob);        

        clawpack46_vt->fort_qinit     = &CLAWPACK46_QINIT;
        clawpack46_vt->fort_setaux    = &CLAWPACK46_SETAUX;
        clawpack46_vt->fort_rpn2      = &CLAWPACK46_RPN2ADV;
        clawpack46_vt->fort_rpt2      = &CLAWPACK46_RPT2ADV;
    }
    else if (user->claw_version == 5)
    {
        fc2d_clawpack5_vtable_t *clawpack5_vt = fc2d_clawpack5_vt(glob);

        clawpack5_vt->fort_qinit     = &CLAWPACK5_QINIT;
        clawpack5_vt->fort_setaux    = &CLAWPACK5_SETAUX;
        clawpack5_vt->fort_rpn2      = &CLAWPACK5_RPN2ADV;
        clawpack5_vt->fort_rpt2      = &CLAWPACK5_RPT2ADV;
    }

    /* Move gauges */
    fclaw_gauges_vtable_t* gauges_vt = (fclaw_gauges_vtable_t*) 
                        fclaw_global_get_vtable(glob, "fclaw_gauges");

    gauges_vt->move = periodic_gauge_move;
}






