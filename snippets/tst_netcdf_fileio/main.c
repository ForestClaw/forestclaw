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

char* concat(const char *s1, const char *s2)
{
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = malloc(len1+len2+1);//+1 for the zero-terminator
    //in real code you would check for errors in malloc here
    memcpy(result, s1, len1);
    memcpy(result+len1, s2, len2+1);//+1 to copy the null-terminator
    return result;
}

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
    double xlower = 0.0;
    double xupper = 1.0;

	const char filename[] = "claw0001.nc";
    const char* dim_name[1] = {"x"};
    char dimnametemp[LEN];
    char* temp;

	int ncid, x_dimid, meqn_dimid, qid;
    int patch0id;
	int dimids[numdim+1];

	int retval;

	/* Create the file. */
   	if ((retval = nc_create(filename, NC_NETCDF4, &ncid)))
    	ERR(retval);

    /* Create subgroup */
    if ((retval = nc_def_grp(ncid, "patch0", &patch0id)))
        ERR(retval);
    // patch0id = ncid;

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
    if ((retval = nc_put_att_string(patch0id, NC_GLOBAL, "dim_names", 1, (const char**) dim_name)))
        ERR(retval);

    // Dimension attribute
    for (int i = 0; i < numdim; ++i)
    {
        strcpy(dimnametemp, dim_name[i]);
        // num_cells
        temp = concat(dimnametemp,".num_cells");
        if ((retval = nc_put_att_int(patch0id, NC_GLOBAL, temp, NC_INT, 1, &mx)))
            ERR(retval);
        // lower bound of this dimension        
        temp = concat(dimnametemp,".lower");
        if ((retval = nc_put_att_double(patch0id, NC_GLOBAL, temp, NC_DOUBLE, 1, &xlower)))
            ERR(retval);
        // upper bound of this dimension
        temp = concat(dimnametemp,".upper");
        if ((retval = nc_put_att_double(patch0id, NC_GLOBAL, temp, NC_DOUBLE, 1, &xupper)))
            ERR(retval);
    }
    if ((retval = nc_put_att_int(patch0id, NC_GLOBAL, "num_eqn", NC_INT, 1, &meqn)))
        ERR(retval);

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
