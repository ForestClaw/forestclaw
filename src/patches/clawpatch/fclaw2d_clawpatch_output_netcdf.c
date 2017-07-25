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

#include <fclaw2d_clawpatch_output_netcdf.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_options.h>

#include <netcdf.h>
#include <netcdf_par.h>

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


#if 0
void cb_clawpatch_output_netcdf_defpatch (fclaw2d_domain_t * domain,
                                          fclaw2d_patch_t * this_patch,
                                          int this_block_idx, int this_patch_idx,
                                          void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;
    int ncid = *((int *) s->user);
    int subgroupid;
    int retval;
    int x_dimid, y_dimid, meqn_dimid, qid;
    int dimids[3];

    int patch_num;
    int level;
    int mx,my,mbc,meqn;
    double xlower,ylower,xupper,yupper,dx,dy;
    double *q;

    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &patch_num,&level);
    
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);


    char patchname[10];
    sprintf(patchname, "patch%04d", patch_num);
    /* Create subgroup */
    if ((retval = nc_def_grp(ncid, patchname, &subgroupid)))
    ERR(retval);

    const char* dim_name[2] = {"x","y"}; 

    xupper = xlower + mx*dx;
    yupper = ylower + my*dy;

    /* Define the dimensions. NetCDF will hand back an ID for each. */
    // dim(q) = mx*my*meqn
    if ((retval = nc_def_dim(subgroupid, "x", mx, &x_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "y", my, &y_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "meqn", meqn, &meqn_dimid)))
        ERR(retval);

    /* Define the attribute */
    // General patch properties
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, "patch_index", NC_INT, 1, &this_patch_idx)))
        ERR(retval);
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, "level", NC_INT, 1, &level)))
        ERR(retval);

    // Dimension name
    if ((retval = nc_put_att_string(subgroupid, NC_GLOBAL, "dim_names", 2, (const char**) dim_name)))
        ERR(retval);

    // Dimension attribute
    // x
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL,"x.num_cells", NC_INT, 1, &mx)))
        ERR(retval);
    if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL,"x.lower", NC_DOUBLE, 1, &xlower)))
        ERR(retval);
    if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL,"x.upper", NC_DOUBLE, 1, &xupper)))
        ERR(retval);
    // y
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL,"y.num_cells", NC_INT, 1, &my)))
        ERR(retval);
    if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL,"y.lower", NC_DOUBLE, 1, &ylower)))
        ERR(retval);
    if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL,"y.upper", NC_DOUBLE, 1, &yupper)))
        ERR(retval);

    dimids[0] = x_dimid;
    dimids[1] = y_dimid;
    dimids[2] = meqn_dimid;

    /* Define the variable. */
    if ((retval = nc_def_var(subgroupid, "q", NC_DOUBLE, 3, dimids, &qid)))
        ERR(retval);

    /* End define mode. */
    if ((retval = nc_enddef(subgroupid)))
        ERR(retval);
}
#endif

#if 0
void cb_clawpatch_output_netcdf_writeq (fclaw2d_domain_t * domain,
                                        fclaw2d_patch_t * this_patch,
                                        int this_block_idx, int this_patch_idx,
                                        void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;

    int ncid = *((int *) s->user);
    int groupid,qid,patchlocationid,patchinfoid;
    int retval;
    int patch_num, level, meqn, mx, my, mbc;
    double xlower, ylower, dx, dy;
    double* q;
    double patchlocation[4];
    int patchinfo[2];

    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &patch_num,&level);
    
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    
    char patchname[10];
    sprintf(patchname, "patch%04d", patch_num);
    if ((retval = nc_inq_grp_ncid(ncid, patchname, &groupid)))
        ERR(retval);
    if ((retval = nc_inq_varid(groupid, "q", &qid)))
        ERR(retval);
    if ((retval = nc_inq_varid(groupid, "patch_location", &patchlocationid)))
        ERR(retval);
    if ((retval = nc_inq_varid(groupid, "patch_info", &patchinfoid)))
        ERR(retval);

    patchinfo[0] = patch_num;
    patchinfo[0] = level;
    if ((retval = nc_put_var_int(groupid, patchinfoid, &patchinfo[0])))
        ERR(retval);

    patchlocation[0] = xlower;
    patchlocation[1] = ylower;
    patchlocation[2] = xlower + dx*mx;
    patchlocation[3] = xlower + dy*my;
    if ((retval = nc_put_var_double(groupid, patchlocationid, &patchlocation[0])))
        ERR(retval);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    if ((retval = nc_put_var_double(groupid, qid, &q[0])))
        ERR(retval);

}

static
void fclaw2d_clawpatch_defgroup_netcdf(fclaw2d_global_t* glob, int ipatch, int ncid)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    
    int subgroupid;
    int retval;
    int numdim = 2;
    int x_dimid, y_dimid, meqn_dimid, location_dimid, info_dimid,
        qid, patchlocationid, patchinfoid;
    int dimids[3];

    int mx,my,meqn;
    meqn = clawpatch_opt->meqn;
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;

    char patchname[10];
    sprintf(patchname, "patch%04d", ipatch);
    /* Create subgroup */
    if ((retval = nc_def_grp(ncid, patchname, &subgroupid)))
    ERR(retval);

    /* Define the dimensions. NetCDF will hand back an ID for each. */
    // dim(q) = mx*my*meqn
    if ((retval = nc_def_dim(subgroupid, "mx", mx, &x_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "my", my, &y_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "meqn", meqn, &meqn_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "mlocation", numdim*2, &location_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "minfo", 2, &info_dimid)))
        ERR(retval);

    dimids[0] = x_dimid;
    dimids[1] = y_dimid;
    dimids[2] = meqn_dimid;

    /* Define the variable. */
    if ((retval = nc_def_var(subgroupid, "patch_info", NC_INT, 1, &info_dimid, &patchinfoid)))
        ERR(retval);

    if ((retval = nc_def_var(subgroupid, "patch_location", NC_DOUBLE, 1, &location_dimid, &patchlocationid)))
        ERR(retval);

    if ((retval = nc_def_var(subgroupid, "q", NC_DOUBLE, 3, dimids, &qid)))
        ERR(retval);

    /* End define mode. */
    if ((retval = nc_enddef(subgroupid)))
        ERR(retval);
}
#endif

void cb_clawpatch_output_netcdf_writeq (fclaw2d_domain_t * domain,
                                        fclaw2d_patch_t * this_patch,
                                        int this_block_idx, int this_patch_idx,
                                        void *user)
{
    fclaw2d_global_iterate_t* s = (fclaw2d_global_iterate_t*) user;
    fclaw2d_global_t *glob = (fclaw2d_global_t*) s->glob;

    int ncid = *((int *) s->user);
    int qid,patchlocationid,patchinfoid;
    int retval;
    size_t start, count;
    int patch_num, level, meqn, mx, my, mbc;
    double xlower, ylower, dx, dy;
    double* q;
    double patchlocation[4];

    /* Get info not readily available to user */
    fclaw2d_patch_get_info(glob->domain,this_patch,
                           this_block_idx,this_patch_idx,
                           &patch_num,&level);
    
    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    
    if ((retval = nc_inq_varid(ncid, "q", &qid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "patch_location", &patchlocationid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, "patch_level", &patchinfoid)))
        ERR(retval);

    if ((retval = nc_var_par_access(ncid, qid, NC_COLLECTIVE)))
        ERR(retval);
    if ((retval = nc_var_par_access(ncid, patchlocationid, NC_COLLECTIVE)))
        ERR(retval);
    if ((retval = nc_var_par_access(ncid, patchinfoid, NC_COLLECTIVE)))
        ERR(retval);

    start = patch_num;
    count = 1;
    if ((retval = nc_put_vara_int(ncid, patchinfoid, &start, &count, &level)))
        ERR(retval);

    patchlocation[0] = xlower;
    patchlocation[1] = ylower;
    patchlocation[2] = xlower + dx*mx;
    patchlocation[3] = xlower + dy*my;
    start = patch_num*4;
    count = 4;
    if ((retval = nc_put_vara_double(ncid, patchlocationid, &start, &count, &patchlocation[0])))
        ERR(retval);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    
    /* Following code only works for clawpack5 */
    for (int m = 0; m < meqn; ++m)
    {
        for (int j = mbc; j < my+mbc; ++j)
        {
            start = patch_num*mx*my*meqn + m*mx*my + (j-mbc)*mx;
            count = mx;
            if ((retval = nc_put_vara_double(ncid, qid, 
                                             &start, 
                                             &count, 
                                             &q[m*(mx+2)*(my+2) + j*(mx+2*mbc) + mbc])))
                ERR(retval);
        }
    }

}
/* This function isn't virtualized;  should it be? */
static
void fclaw2d_clawpatch_header_netcdf(fclaw2d_global_t* glob, int iframe, int ncid)
{

    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    int numdim = 2;
    int retval;
    int meqn,ngrids,maux,mx,my;

    double time = glob->curr_time;

    ngrids = glob->domain->global_num_patches;

    meqn = clawpatch_opt->meqn;
    maux = clawpatch_opt->maux;
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;


    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "ngrids", NC_INT, 1, &ngrids)))
        ERR(retval);
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "dim", NC_INT, 1, &numdim)))
        ERR(retval);
    if ((retval = nc_put_att_double(ncid, NC_GLOBAL, "t", NC_DOUBLE, 1, &time)))
        ERR(retval);
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "meqn", NC_INT, 1, &meqn)))
        ERR(retval);
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "maux", NC_INT, 1, &maux)))
        ERR(retval);
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "mx", NC_INT, 1, &mx)))
        ERR(retval);
    if ((retval = nc_put_att_int(ncid, NC_GLOBAL, "my", NC_INT, 1, &my)))
        ERR(retval);    
}

static
void fclaw2d_clawpatch_defgroup_netcdf(fclaw2d_global_t* glob, int ngrids, int ncid)
{
    const fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);
    
    int retval;
    int numdim = 2;
    int q_dimid, location_dimid, level_dimid,
        qid, patchlocationid, patchlevelid;

    int mx,my,meqn;
    meqn = clawpatch_opt->meqn;
    mx = clawpatch_opt->mx;
    my = clawpatch_opt->my;


    /* Define the dimensions. NetCDF will hand back an ID for each. */
    // dim(q) = mx*my*meqn
    if ((retval = nc_def_dim(ncid, "ngridsxmxxmyxmeqn", mx*my*meqn*ngrids, &q_dimid)))
        ERR(retval);
    
    if ((retval = nc_def_dim(ncid, "ngridsxmlocation", numdim*2*ngrids, &location_dimid)))
        ERR(retval);
    
    if ((retval = nc_def_dim(ncid, "ngrids", ngrids, &level_dimid)))
        ERR(retval);

    /* Define the variable. */
    if ((retval = nc_def_var(ncid, "patch_level", NC_INT, 1, &level_dimid, &patchlevelid)))
        ERR(retval);

    if ((retval = nc_def_var(ncid, "patch_location", NC_DOUBLE, 1, &location_dimid, &patchlocationid)))
        ERR(retval);

    if ((retval = nc_def_var(ncid, "q", NC_DOUBLE, 1, &q_dimid, &qid)))
        ERR(retval);

    /* End define mode. */
    if ((retval = nc_enddef(ncid)))
        ERR(retval);
}
    /*--------------------------------------------------------------------
    Public interface
    Use this function as follows : 
           fclaw2d_vtable_t *vt = fclaw2d_vt();
           vt->output_frame = &fclaw2d_clawpatch_output_netcdf;
    -------------------------------------------------------------------- */

void fclaw2d_clawpatch_output_netcdf(fclaw2d_global_t* glob,int iframe)
{
    fclaw2d_domain_t *domain = glob->domain;
    int retval;
    int ncid;
    int ngrids;

    /* BEGIN NON-SCALABLE CODE */
    /* Write the file contents in serial.
       Use only for small numbers of processors. */
    // fclaw2d_domain_serialization_enter (domain);

    /* Create the file. */
    char filename[12];
    sprintf(filename,"claw%04d.nc",iframe);

    /* serial */
    // if ((retval = nc_create(filename, NC_NETCDF4, &ncid)))
        // ERR(retval);

    /* parallel */
    if ((retval = nc_create_par(filename, NC_NETCDF4|NC_MPIIO, glob->mpicomm, 
                                MPI_INFO_NULL, &ncid)));

    fclaw2d_clawpatch_header_netcdf(glob, iframe, ncid);
    
    /* Create groups for patches */
    ngrids = domain->global_num_patches;
    fclaw2d_clawpatch_defgroup_netcdf(glob,ngrids,ncid);

#if 0
    for (int i = 0; i < ngrids; ++i)
    {
        fclaw2d_clawpatch_defgroup_netcdf(glob,i,ncid);
    }
#endif

    /* Write out patch data*/
    fclaw2d_global_iterate_patches (glob, cb_clawpatch_output_netcdf_writeq, &ncid);
    /* Close the file. */
    if ((retval = nc_close(ncid)))
        ERR(retval);
    // fclaw2d_domain_serialization_leave (domain);
    /* END OF NON-SCALABLE CODE */
}

