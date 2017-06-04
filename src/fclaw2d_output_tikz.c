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

#include <fclaw2d_output.h>
#include <fclaw_base.h>   /* Needed for MPI declarations */

#include <fclaw2d_domain.h>
#include <fclaw2d_global.h>

#include <fclaw2d_options.h>
#include <fclaw_math.h>

typedef struct
{
    int mx;
    int my;
    double ax;
    double bx;
    double ay;
    double by;
    FILE *fp;

} fclaw2d_tikz_info_t;

static void
cb_tikz_output (fclaw2d_domain_t * domain,
                fclaw2d_patch_t * this_patch,
                int this_block_idx, int this_patch_idx,
                void *user)
{
    fclaw2d_global_iterate_t *s = (fclaw2d_global_iterate_t *) user;
    fclaw2d_tikz_info_t *s_tikz = (fclaw2d_tikz_info_t*) s->user;

    int mx, my, mi, mj, level, lmax, mxf, myf;
    double ax,bx,ay,by,dxf,dyf;
    int xlow_d, ylow_d, xupper_d, yupper_d;

    FILE *fp = s_tikz->fp;
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(s->glob);

    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num = domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);

    level = this_patch->level;
    lmax = fclaw_opt->maxlevel;

    mi = fclaw_opt->mi;
    mj = fclaw_opt->mj;

    mx = s_tikz->mx;
    my = s_tikz->my;

    ax = s_tikz->ax;
    ay = s_tikz->ay;
    bx = s_tikz->bx;
    by = s_tikz->by;

    mxf = mi*mx*pow_int(2,lmax-level);
    myf = mj*my*pow_int(2,lmax-level);

    dxf = 1.0/(mi*mx*pow_int(2,lmax));
    dyf = 1.0/(mj*my*pow_int(2,lmax));

    xlow_d = (this_patch->xlower-ax)/dxf;
    ylow_d = (this_patch->ylower-ay)/dyf;
    xupper_d = xlow_d + mxf;
    yupper_d = ylow_d + myf;

    fprintf(fp,"    %% Patch number %ld; rank = %d\n",(long int) patch_num,
            s->glob->domain->mpirank);
    fprintf(fp,"    \\draw [ultra thin] (%d,%d) rectangle (%d,%d);\n\n",
            xlow_d,ylow_d,xupper_d,yupper_d);
}

void fclaw2d_output_frame_tikz(fclaw2d_global_t* glob, int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;
    fclaw2d_tikz_info_t s_tikz;

    char fname[20];

    /* Should be in fclaw_opt */
    const fclaw_options_t *fclaw_opt = fclaw2d_get_options(glob);
    
    double figsize[2];
    figsize[0] = fclaw_opt->tikz_figsize[0];   /* Inches */
    figsize[1] = fclaw_opt->tikz_figsize[1];   /* Inches */

    s_tikz.mx = 1;  /* Doesn't really matter what we put here;  used only for scaling */
    s_tikz.my = 1;

    if (!fclaw_opt->manifold)
    {
        s_tikz.ax = fclaw_opt->ax;
        s_tikz.bx = fclaw_opt->bx;
        s_tikz.ay = fclaw_opt->ay;
        s_tikz.by = fclaw_opt->by;
    }
    else
    {
        /* Can't do much with a manifold */
        s_tikz.ax = 0;
        s_tikz.bx = 1;
        s_tikz.ay = 0;
        s_tikz.by = 1;
    }

    int lmax = fclaw_opt->maxlevel;
    int mi = fclaw_opt->mi;
    int mj = fclaw_opt->mj;

    int mxf = mi*s_tikz.mx*pow_int(2,lmax);  
    int myf = mj*s_tikz.my*pow_int(2,lmax);
    double sx = figsize[0]/mxf;
    double sy = figsize[1]/myf;

    FILE *fp;
    sprintf(fname,"tikz_%04d.tex",iframe);  /* fname[20] */
    if (domain->mpirank == 0)
    {
        /* Only rank 0 opens the file */
        fp = fopen(fname,"w");

        fprintf(fp,"\\documentclass{standalone}\n");
        fprintf(fp,"\\usepackage{tikz}\n");
        fprintf(fp,"\n");
        fprintf(fp,"\\newcommand{\\plotfig}[1]{#1}\n");
        fprintf(fp,"\\newcommand{\\plotgrid}[1]{#1}\n");
        fprintf(fp,"\\newcommand{\\figname}{%s_%04d.%s}\n",fclaw_opt->tikz_plot_prefix,
                iframe,fclaw_opt->tikz_plot_suffix);
        fprintf(fp,"\n");
        fprintf(fp,"\\begin{document}\n");
        fprintf(fp,"\\begin{tikzpicture}[x=%18.16fin, y=%18.16fin]\n",sx,sy);
        fprintf(fp,"    \\plotfig{\\node (forestclaw_plot) at (%3.1f,%3.1f)\n",
                ((double) mxf)/2,((double) myf)/2);
        fprintf(fp,"    {\\includegraphics{\\figname}};}\n\n");
        fprintf(fp,"\\plotgrid{\n");
        fclose(fp); 
    }

    fp = fopen(fname,"a"); 
    s_tikz.fp = fp;

    fclaw2d_domain_serialization_enter (domain);
    fclaw2d_global_iterate_patches (glob, cb_tikz_output, (void *) &s_tikz);
    fclaw2d_domain_serialization_leave (domain);

    fclose(fp);

    /* Wait for all processes to finish before adding last lines of file */
    
#ifdef FCLAW_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (domain->mpirank == 0)
    {
        fp = fopen(fname,"a");
        fprintf(fp,"} %% end plotgrid\n");
        fprintf(fp,"\\end{tikzpicture}\n");
        fprintf(fp,"\\end{document}\n");
        fclose(fp);
    }
}
