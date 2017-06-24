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
	int num_aux = 0;
	int num_dim = 1;
	int meqn = 1;
	int mx = 5;
	double q[] = {0.1,0.2,0.3,0.4,0.5};

	const char filename[] = "claw0000.nc";
	int ncid, x_dimid, meqn_dimid, qid;
	int dimids[num_dim+1];

	int retval;

	/* Create the file. */
   	if ((retval = nc_create(filename, NC_CLOBBER, &ncid)))
    	ERR(retval);
    
    /* Define the dimensions. NetCDF will hand back an ID for each. */
    if ((retval = nc_def_dim(ncid, "x", mx, &x_dimid)))
    	ERR(retval);
    if ((retval = nc_def_dim(ncid, "meqn", meqn, &meqn_dimid)))
    	ERR(retval);

    dimids[0] = x_dimid;
    dimids[1] = meqn_dimid;

    /* Define the variable. */
    if ((retval = nc_def_var(ncid, "q", NC_DOUBLE, num_dim+1, dimids, &qid)))
    	ERR(retval);
    
    /* End define mode. */
    if ((retval = nc_enddef(ncid)))
    	ERR(retval);

    /* Write the pretend data to the file. Although netCDF supports
     * reading and writing subsets of data, in this case we write all
     * the data in one operation. */
    if ((retval = nc_put_var_double(ncid, qid, &q[0])))
    	ERR(retval);

    /* Close the file. */
    if ((retval = nc_close(ncid)))
    	ERR(retval);

    printf("*** SUCCESS writing file claw0000.nc!\n");

	return 0;
}
