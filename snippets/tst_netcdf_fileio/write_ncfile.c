//
// This file is modified from netcdf-4.4.1.1/examples/simple_xy_wr.c
//


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define LEN 40

void write_patch(int ncid, int subgroupid, int patchidx,
                 double time, int ngrids, int numdim, int meqn, int maux, int mx)
{
    const char* dim_name[1] = {"x"};
    // char* temp;
    char temp[LEN];

    /* Patch specific */
    int level = 0;
    double q[mx];
    double xlower = 0.0;
    double xupper = 1.0;

    int x_dimid, meqn_dimid, qid;
    int dimids[numdim+1];

    int retval;

    /* Create pretend data*/
    for (int i = 0; i < mx; ++i)
    {
        q[i] = patchidx + i;
    }

    /* Define the dimensions. NetCDF will hand back an ID for each. */
    // dim(q) = mx*my*meqn 
    if ((retval = nc_def_dim(subgroupid, "x", mx, &x_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(subgroupid, "meqn", meqn, &meqn_dimid)))
        ERR(retval);

    /* Define the attribute */
    // General patch properties
    if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL, "t", NC_DOUBLE, 1, &time)))
        ERR(retval);
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, "num_eqn", NC_INT, 1, &meqn)))
        ERR(retval);
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, "num_aux", NC_INT, 1, &maux)))
        ERR(retval);
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, "patch_index", NC_INT, 1, &patchidx)))
        ERR(retval);
    if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, "level", NC_INT, 1, &level)))
        ERR(retval);
    
    // Dimension name
    if ((retval = nc_put_att_string(subgroupid, NC_GLOBAL, "dim_names", 1, (const char**) dim_name)))
        ERR(retval);

    // Dimension attribute
    for (int i = 0; i < numdim; ++i)
    {
        // num_cells
        sprintf(temp,"%s.num_cells",dim_name[i]);
        if ((retval = nc_put_att_int(subgroupid, NC_GLOBAL, temp, NC_INT, 1, &mx)))
            ERR(retval);

        // lower bound of this dimension        
        sprintf(temp,"%s.lower",dim_name[i]);
        if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL, temp, NC_DOUBLE, 1, &xlower)))
            ERR(retval);

        // upper bound of this dimension
        sprintf(temp,"%s.upper",dim_name[i]);
        if ((retval = nc_put_att_double(subgroupid, NC_GLOBAL, temp, NC_DOUBLE, 1, &xupper)))
            ERR(retval);
    }

    dimids[0] = x_dimid;
    dimids[1] = meqn_dimid;
    /* Define the variable. */
    if ((retval = nc_def_var(subgroupid, "q", NC_DOUBLE, numdim+1, dimids, &qid)))
        ERR(retval);
    
    /* End define mode. */
    if ((retval = nc_enddef(subgroupid)))
        ERR(retval);

    /* Write the pretend data to the file. Although netCDF supports
     * reading and writing subsets of data, in this case we write all
     * the data in one operation. */
    if ((retval = nc_put_var_double(subgroupid, qid, &q[0])))
        ERR(retval);
}

int main(int argc, char const *argv[])
{
    const char filename[] = "claw0001.nc";
    int ncid;
    int numpatch = 2;
    int subgroupid[numpatch];
    int retval;
    char patchname[7];

    /* Variables that general to all patches */
    double time = 0.0;
    int ngrids = 1;
    int numdim = 1;
    int meqn = 1;
    int maux = 0;
    int mx = 5;

    /* Create the file. */
    if ((retval = nc_create(filename, NC_NETCDF4, &ncid)))
        ERR(retval);

    for (int i = 0; i < numpatch; ++i)
    {
        sprintf(patchname, "patch%d", i);
        /* Create subgroup */
        if ((retval = nc_def_grp(ncid, patchname, &subgroupid[i])))
            ERR(retval);

        write_patch(ncid, subgroupid[i], i, time, ngrids, numdim, meqn, maux, mx);
    }
    /* Close the file. */
    if ((retval = nc_close(ncid)))
        ERR(retval);

    printf("*** SUCCESS writing file %s!\n", filename);

    return 0;
}