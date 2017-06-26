//
// This file is modified from netcdf-4.4.1.1/examples/simple_xy_wr.c
//

#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main(int argc, char const *argv[])
{
	/* state variable */
	double time = 0.0;
	int ngrids = 1;
	int maux = 0;
	int numdim = 1;
	int meqn = 1;
	int mx = 5;
    int patchidx = 1;
    int level = 0;
	double q[] = {0.1,0.2,0.3,0.4,0.5};

	const char filename[] = "claw0001.nc";
    const char* dim_name[1] = {"x"};
    const char* greeting[2] = {"Hello\0", "world\0"};

	int ncid, x_dimid, meqn_dimid, qid;
    int patch0id;
	int dimids[numdim+1];

	int retval;

	/* Create the file. */
   	if ((retval = nc_create(filename, NC_CLOBBER, &ncid)))
    	ERR(retval);

    /* Create subgroup */
    // if ((retval = nc_def_grp(ncid, "patch0", &patch0id)))
    //     ERR(retval);
    patch0id = ncid;

    /* Define the attribute */
    // General patch properties
    if ((retval = nc_put_att_double(patch0id, NC_GLOBAL, "t", NC_DOUBLE, 1, &time)))
        ERR(retval);
    if ((retval = nc_put_att_int(patch0id, NC_GLOBAL, "num_eqn", NC_INT, 1, &meqn)))
        ERR(retval);
    if ((retval = nc_put_att_int(patch0id, NC_GLOBAL, "num_aux", NC_INT, 1, &maux)))
        ERR(retval);
    if ((retval = nc_put_att_int(patch0id, NC_GLOBAL, "patch_index", NC_INT, 1, &patchidx)))
        ERR(retval);
    if ((retval = nc_put_att_int(patch0id, NC_GLOBAL, "level", NC_INT, 1, &level)))
        ERR(retval);
    
    // Dimension name
    // Following code is not yet working ...
    // if ((retval = nc_put_att_string(patch0id, NC_GLOBAL, "dim_names", 1, (const char**) dim_name)))
    //     printf("%d\n", retval);
    //     ERR(retval);

    // Dimension attribute

    /* Define the dimensions. NetCDF will hand back an ID for each. */
    // dim(q) = mx*my*meqn 
    if ((retval = nc_def_dim(patch0id, "x", mx, &x_dimid)))
    	ERR(retval);
    if ((retval = nc_def_dim(patch0id, "meqn", meqn, &meqn_dimid)))
    	ERR(retval);

    dimids[0] = x_dimid;
    dimids[1] = meqn_dimid;

    /* Define the variable. */
    if ((retval = nc_def_var(patch0id, "q", NC_DOUBLE, numdim+1, dimids, &qid)))
    	ERR(retval);
    
    /* End define mode. */
    if ((retval = nc_enddef(patch0id)))
    	ERR(retval);

    /* Write the pretend data to the file. Although netCDF supports
     * reading and writing subsets of data, in this case we write all
     * the data in one operation. */
    if ((retval = nc_put_var_double(patch0id, qid, &q[0])))
    	ERR(retval);

    /* Close the file. */
    if ((retval = nc_close(ncid)))
    	ERR(retval);

    printf("*** SUCCESS writing file %s!\n", filename);

	return 0;
}
