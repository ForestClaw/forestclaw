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

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_output.h>
#include <fclaw2d_vtk.h>
#include <fclaw2d_map.h>

#include <fclaw2d_clawpatch.h>


static void
cb_serial_output (fclaw2d_domain_t * domain,
                  fclaw2d_patch_t * this_patch,
                  int this_block_idx, int this_patch_idx,
                  void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;

    int iframe = *((int *) s->user);
    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
        domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);

    /* TODO Enable 64bit integers for global counters and indices */

    /* the user can also get this, but maybe we don't want the user
       to have access? */
    int level = this_patch->level;

    fclaw2d_patch_write_file(s->glob, this_patch, this_block_idx,
                             this_patch_idx, iframe, (int) patch_num,
                             level);
}

static
void fclaw2d_output_write_serial(fclaw2d_global_t* glob,int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    fclaw2d_domain_serialization_enter (domain);

    if (domain->mpirank == 0)
    {
        fclaw2d_patch_write_header(glob,iframe);
    }

    /* Write out each patch to fort.qXXXX */
    fclaw2d_global_iterate_patches (glob, cb_serial_output, (void *) &iframe);
    fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}


static void
cb_tikz_output (fclaw2d_domain_t * domain,
                fclaw2d_patch_t * this_patch,
                int this_block_idx, int this_patch_idx,
                void *user)
{
    fclaw2d_global_iterate_t *s = (fclaw2d_global_iterate_t *) user;

    FILE *fp = (FILE*) s->user;
    const amr_options_t *gparms = s->glob->gparms;

    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    int64_t patch_num =
        domain->global_num_patches_before +
        (int64_t) (this_block->num_patches_before + this_patch_idx);

    int level = this_patch->level;
    int lmax = gparms->maxlevel;

    int mx, my, mbc;
    double dx,dy;
    double xlower, ylower;
    fclaw2d_clawpatch_grid_data(s->glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double ax,bx,ay,by,dxf,dyf;
    if (!gparms->manifold)
    {
        ax = gparms->ax;
        bx = gparms->bx;
        ay = gparms->ay;
        by = gparms->by;
        int mi = gparms->mi;
        int mj = gparms->mj;
        dxf = (bx-ax)/(mi*mx*pow_int(2,lmax));
        dyf = (by-ay)/(mj*my*pow_int(2,lmax));
    }
    else
    {
        ax = 0;
        ay = 0;
        dxf = 1.0/(mx*pow_int(2,lmax));
        dyf = 1.0/(my*pow_int(2,lmax));
    }


    int mxf = mx*pow_int(2,lmax-level);
    int myf = my*pow_int(2,lmax-level);

    int xlow_d = (xlower-ax)/dxf;
    int ylow_d = (ylower-ay)/dyf;
    int xupper_d = xlow_d + mxf;
    int yupper_d = ylow_d + myf;

    fprintf(fp,"    %% Patch number %ld\n",(long int) patch_num);
    fprintf(fp,"    \\draw [ultra thin] (%d,%d) rectangle (%d,%d);\n\n",
            xlow_d,ylow_d,xupper_d,yupper_d);
}


static
void fclaw2d_output_write_tikz(fclaw2d_global_t* glob,int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;

    char fname[20];
    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */

    /* Should be in gparms */
    const amr_options_t *gparms = glob->gparms;
    double figsize[2];
    figsize[0] = gparms->tikz_figsize[0];   /* Inches */
    figsize[1] = gparms->tikz_figsize[1];   /* Inches */

    int lmax = gparms->maxlevel;
    int mx = gparms->mx;
    int my = gparms->my;
    int mi = gparms->mi;
    int mj = gparms->mj;

    int mxf = mi*mx*pow_int(2,lmax);
    int myf = mj*my*pow_int(2,lmax);
    double sx = figsize[0]/mxf;
    double sy = figsize[1]/myf;

#ifdef FCLAW_ENABLE_MPI

    /* Don't do anything until we are done with this part */
    MPI_Barrier(MPI_COMM_WORLD);
    FILE *fp;
    sprintf(fname,"tikz.%04d.tex",iframe);  /* fname[20] */
    if (domain->mpirank == 0)
    {
        /* Only rank 0 opens the file */
        fp = fopen(fname,"w");
        fprintf(fp,"\\begin{tikzpicture}[x=%18.16fin, y=%18.16fin]\n",sx,sy);
        fprintf(fp,"    \\node (forestclaw_plot) at (%3.1f,%3.1f)\n",
                ((double) mxf)/2,((double)myf)/2);
        fprintf(fp,"    {\\includegraphics{\\figname}};\n\n");
        fprintf(fp,"\\plotgrid{\n");
        fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    int i;
    for(i = 0; i < domain->mpisize; i++)
    {
        if (domain->mpirank == i)
        {
            fp = fopen(fname,"a");
            fclaw2d_global_iterate_patches (glob, cb_tikz_output, (void *) fp);
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (domain->mpirank == 0)
    {
        fp = fopen(fname,"a");
        fprintf(fp,"} %% end plotgrid\n");
        fprintf(fp,"\\end{tikzpicture}\n");
        fclose(fp);
    }

    /* This might not be needed ... */
    MPI_Barrier(MPI_COMM_WORLD);

    /* END OF NON-SCALABLE CODE */
#else
#if 0
    FILE *fp;
    sprintf(fname,"tikz.%04d.tex",iframe);  /* fname[20] */
    /* Only rank 0 opens the file */
    fp = fopen(fname,"w");
    fprintf(fp,"\\begin{tikzpicture}[x=%18.16fin, y=%18.16fin]\n",sx,sy);
    fprintf(fp,"    \\node (forestclaw_plot) at (%3.1f,%3.1f)\n",double(mxf)/2,double(myf)/2);
    fprintf(fp,"    {\\includegraphics{\\figname}};\n\n");
    fprintf(fp,"\\plotgrid{\n");

    /* Write out each patch to tikz file */
    fclaw2d_domain_iterate_patches (domain, cb_tikz_output, (void *) fp);

    fp = fopen(fname,"a");
    fprintf(fp,"} %% end plotgrid\n");
    fprintf(fp,"\\end{tikzpicture}\n");
    fclose(fp);
#endif
#endif
}


static void
fclaw2d_output_vtk_coordinate_cb (fclaw2d_global_t * glob,
                                  fclaw2d_patch_t * this_patch,
                                  int this_block_idx, int this_patch_idx,
                                  char *a)
{
    const amr_options_t *gparms = glob->gparms;
    fclaw2d_map_context_t *cont;

    int mx,my,mbc;
    double dx,dy,xlower,ylower;


    cont = glob->cont;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

#if 0
    fclaw2d_clawpatch_metric_data(glob,this_patch,&xp,&yp,&zp,&xd,&yd,&zd,&area);

    ClawPatch *cp = fclaw2d_clawpatch_get_cp(this_patch);
    const double xlower = cp->xlower ();
    const double ylower = cp->ylower ();
    const double dx = cp->dx ();
    const double dy = cp->dy ();
#endif

    /* Enumerate point coordinates in the patch */
    double *d = (double *) a;
    int i, j;
    double xpp,ypp,zpp;
    for (j = 0; j <= my; ++j)
    {
        const double y = ylower + j * dy;
        for (i = 0; i <= mx; ++i)
        {
            const double x = xlower + i * dx;
            if (gparms->manifold)
            {
                FCLAW2D_MAP_C2M(&cont,&this_block_idx,&x,&y,&xpp,&ypp,&zpp);
                *d++ = xpp;
                *d++ = ypp;
                *d++ = zpp;
            }
            else
            {
                *d++ = x;
                *d++ = y;
                *d++ = 0;
            }
        }
    }
}


static void
fclaw2d_output_vtk_value_cb (fclaw2d_global_t * glob,
                             fclaw2d_patch_t * this_patch,
                             int this_block_idx, int this_patch_idx,
                             char *a)
{
    double *q;
    double xlower,ylower,dx,dy;
    int mx,my,mbc,meqn;

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    const int xlane = mx + 2 * mbc;
    const int ylane = my + 2 * mbc;

    // Enumerate equation data in the patch
    float *f = (float *) a;
    int i, j, k;
    for (j = 0; j < my; ++j)
    {
        for (i = 0; i < mx; ++i)
        {
            for (k = 0; k < meqn; ++k)
            {
                *f++ = (float) q[(k * ylane + j + mbc) * xlane + i + mbc];
            }
        }
    }
}

void
fclaw2d_output_write_vtk (fclaw2d_global_t * glob, const char *basename)
{
    const amr_options_t *gparms = glob->gparms;

    (void) fclaw2d_vtk_write_file (glob, basename,
                                   gparms->mx, gparms->my, gparms->meqn,
                                   gparms->vtkspace, gparms->vtkwrite,
                                   fclaw2d_output_vtk_coordinate_cb,
                                   fclaw2d_output_vtk_value_cb);
}

void
fclaw2d_output_frame (fclaw2d_global_t * glob, int iframe)
{
    double time;

    time = glob->curr_time;

    /* Record output time */
    fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_OUTPUT]);

    /* Output VTK file while we're at it */

    const amr_options_t *gparms = glob->gparms;
    if (gparms->vtkout & 2)
    {
        if (gparms->serialout)
            fclaw_global_essentialf("Output Frame %4d  at time %16.8e (vtk,ascii)\n\n",
                                    iframe,time);
        else
            fclaw_global_essentialf("Output Frame %4d  at time %16.8e (vtk)\n\n",
                                    iframe,time);
    }
    else if (gparms->serialout)
    {
        fclaw_global_essentialf("Output Frame %d  at time %16.8e (ascii)\n\n",
                                iframe,time);
    }
    else
    {
        fclaw_global_essentialf("Output Frame %d  at time %16.8e (none)\n\n",
                                iframe,time);
    }

    if (gparms->vtkout & 2)
    {
        char basename[BUFSIZ];
        snprintf (basename, BUFSIZ, "%s_frame_%04d", gparms->prefix, iframe);
        fclaw2d_output_write_vtk (glob, basename);
    }

    if (gparms->serialout)
    {
        fclaw2d_output_write_serial(glob, iframe);
    }

    if (gparms->tikzout)
    {
        fclaw2d_timer_start (&glob->timers[FCLAW2D_TIMER_EXTRA3]);
        fclaw2d_output_write_tikz(glob,iframe);
        fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_EXTRA3]);
    }

    /* Record output time */
    fclaw2d_timer_stop (&glob->timers[FCLAW2D_TIMER_OUTPUT]);
}
