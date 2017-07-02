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

void read_patch(int ncid, int subgroupid, int numdim)
{
    char* dim_name[numdim];
    // char* temp;
    char temp[LEN];

    /* Patch specific */
    int level = 0;
    double* q;
    double xlower;
    double xupper;
    double time;
    int meqn, maux, patchidx, mx;

    int qid;
    int dimids[numdim+1];

    int retval;

    /* Define the attribute */
    // General patch properties
    if ((retval = nc_get_att_double(subgroupid, NC_GLOBAL, "t", &time)))
        ERR(retval);
    if ((retval = nc_get_att_int(subgroupid, NC_GLOBAL, "num_eqn", &meqn)))
        ERR(retval);
    if ((retval = nc_get_att_int(subgroupid, NC_GLOBAL, "num_aux", &maux)))
        ERR(retval);
    if ((retval = nc_get_att_int(subgroupid, NC_GLOBAL, "patch_index", &patchidx)))
        ERR(retval);
    if ((retval = nc_get_att_int(subgroupid, NC_GLOBAL, "level", &level)))
        ERR(retval);
    
    // Dimension name
    if ((retval = nc_get_att_string(subgroupid, NC_GLOBAL, "dim_names", dim_name)))
        ERR(retval);

    // Dimension attribute
    for (int i = 0; i < numdim; ++i)
    {
        // num_cells
        sprintf(temp,"%s.num_cells",dim_name[i]);
        if ((retval = nc_get_att_int(subgroupid, NC_GLOBAL, temp, &mx)))
            ERR(retval);

        // lower bound of this dimension        
        sprintf(temp,"%s.lower",dim_name[i]);
        if ((retval = nc_get_att_double(subgroupid, NC_GLOBAL, temp, &xlower)))
            ERR(retval);

        // upper bound of this dimension
        sprintf(temp,"%s.upper",dim_name[i]);
        if ((retval = nc_get_att_double(subgroupid, NC_GLOBAL, temp, &xupper)))
            ERR(retval);
    }
    /* Allocate memory for q */
    q = malloc(sizeof(double)*mx);
    /* Inquire variable id */
    if ((retval = nc_inq_varid(subgroupid, "q", &qid)))
        ERR(retval);

    /* Get the pretend data to the file.*/
    if ((retval = nc_get_var_double(subgroupid, qid, &q[0])))
        ERR(retval);
    
    for (int i = 0; i < mx; ++i)
    {
        printf("%f\n", q[i]);
    }
    printf("\n");
}

int main(int argc, char const *argv[])
{
    const char filename[] = "claw0001.nc";
    int ncid;
    int numpatch;
    int subgroupid;
    int retval;
    char patchname[7];

    /* Variables that general to all patches */
    double time = 0.0;
    int ngrids = 1;
    int numdim = 1;
    int meqn = 1;
    int maux = 0;
    int mx = 5;

    /* Open the file. */
    if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
        ERR(retval);

    if ((retval = nc_get_att_int(ncid, NC_GLOBAL, "num_patch", &numpatch)))
            ERR(retval);
    if ((retval = nc_get_att_int(ncid, NC_GLOBAL, "num_dim", &numdim)))
            ERR(retval);

    for (int i = 0; i < numpatch; ++i)
    {
        sprintf(patchname, "patch%d", i);
        /* Inquire subgroup id*/
        if ((retval = nc_inq_ncid(ncid, patchname, &subgroupid)))
            ERR(retval);
        printf("Reading %s\n", patchname);
        read_patch(ncid, subgroupid, numdim);
    }

    /* Close the file. */
    if ((retval = nc_close(ncid)))
        ERR(retval);

    printf("*** SUCCESS reading file %s!\n", filename);

    return 0;
}