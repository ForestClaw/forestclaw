#include <fclaw_clawpatch_output_vtpd.h>
#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_filesystem.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.h>
#include <fclaw_math.h>
#include <stdio.h>
#include <stdint.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void
WriteVTIHeader(FILE* file, fclaw_global_t * glob, fclaw_patch_t * patch,
     int blockno, int patchno)
{
    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;
    fclaw_clawpatch_3d_grid_data(glob, patch, 
        &mx, &my, &mz, &mbc, 
        &xlower, &ylower, &zlower, 
        &dx, &dy, &dz);

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(file, "\t<ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%f %f %f\" Spacing=\"%f %f %f\" Direction=\"1 0 0 0 1 0 0 0 1\">\n", mx, my, mz, xlower, ylower, zlower, dx, dy, dz);
    fprintf(file, "\t\t<Piece Extent=\"0 %d 0 %d 0 %d\">\n", mx, my, mz);
    fprintf(file, "\t\t\t<CellData>\n");
}

static uint32_t
WriteQDataArray(FILE* file, fclaw_global_t* glob, fclaw_patch_t* patch)
{
    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;
    fclaw_clawpatch_3d_grid_data(glob, patch, 
        &mx, &my, &mz, &mbc, 
        &xlower, &ylower, &zlower, 
        &dx, &dy, &dz);

    int meqn;
    double* q;
    fclaw_clawpatch_soln_data(glob, patch, &q, &meqn);

    double min_val = INFINITY;
    double max_val = -min_val;

    int stride_j = mx + 2 * mbc;
    int stride_k = stride_j * (my + 2 * mbc);
    int stride_m = stride_k * (mz + 2 * mbc);

    for(int m = 0; m < meqn; m++)
    for(int k = 0; k < mz; k++)
    for(int j = 0; j < my; j++)
    for(int i = 0; i < my; i++)
    {
        double val = q[mbc + i + (j+mbc) * stride_j + (k+mbc) * stride_k+m*stride_m];
        min_val = MIN(min_val, val);
        max_val = MAX(max_val, val);
    }


    fprintf(file, "\t\t\t\t<DataArray type=\"Float32\" Name=\"%s\" NumberOfComponents=\"%d\" format=\"appended\" RangeMin=\"%f\" RangeMax=\"%f\" offset=\"0\">\n", 
                          "meqn", meqn, min_val, max_val);
    fprintf(file, "\t\t\t\t</DataArray>\n");

    uint32_t offset = sizeof(uint32_t) + mx * my * mz * meqn * sizeof(float);
    return offset;
}

static void
WriteIntDataArray(FILE* file, int offset, const char* name, int value)
{
    fprintf(file, "\t\t\t\t<DataArray type=\"Int32\" Name=\"%s\" NumberOfComponents=\"1\" format=\"appended\" RangeMin=\"%d\" RangeMax=\"%d\" offset=\"%d\">\n", 
                          name, value, value, offset);
    fprintf(file, "\t\t\t\t</DataArray>\n");
}

static void
AppendQData(FILE* file, fclaw_global_t* glob, fclaw_patch_t* patch)
{
    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;
    fclaw_clawpatch_3d_grid_data(glob, patch, 
        &mx, &my, &mz, &mbc, 
        &xlower, &ylower, &zlower, 
        &dx, &dy, &dz);

    int meqn;
    double* q;
    fclaw_clawpatch_soln_data(glob, patch, &q, &meqn);

    int stride_j = mx + 2 * mbc;
    int stride_k = stride_j * (my + 2 * mbc);
    int stride_m = stride_k * (mz + 2 * mbc);

    uint32_t blocksize = mx * my * mz * meqn * sizeof(float);

    fwrite(&blocksize, sizeof(blocksize), 1, file);

    for(int k = 0; k < mz; k++)
    for(int j = 0; j < my; j++)
    for(int i = 0; i < my; i++)
    for(int m = 0; m < meqn; m++)
    {
        float val = q[i+mbc + (j+mbc) * stride_j + (k+mbc) * stride_k + m * stride_m];
        fwrite(&val, sizeof(float), 1, file);
    }
}

static void
AppendIntData(FILE* file, int value, uint32_t data_length)
{
    fwrite(&data_length, sizeof(data_length), 1, file);
    for(uint32_t i=0; i < data_length/sizeof(int); i++)
    {
        fwrite(&value, sizeof(value), 1, file);
    }
}

static void 
WriteVTIFile(fclaw_domain_t * domain, fclaw_patch_t * patch,
     int blockno, int patchno, void *user)
{
    fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;
    fclaw_global_t* glob = g->glob;
    char *base_name = (char*) g->user;

    const int globalno = glob->domain->global_num_patches_before + glob->domain->blocks[blockno].num_patches_before + patchno;

    char file_name[BUFSIZ];
    snprintf(file_name, BUFSIZ, "%s/%s_%d.vti", base_name, base_name, globalno);

    FILE *file = fopen(file_name, "w");

    int mx, my, mz, mbc;
    double xlower, ylower, zlower, dx, dy, dz;
    fclaw_clawpatch_3d_grid_data(glob, patch, 
        &mx, &my, &mz, &mbc, 
        &xlower, &ylower, &zlower, 
        &dx, &dy, &dz);


    uint32_t int_data_length = mx*my*mz*sizeof(int);

    // Write header
    WriteVTIHeader(file, glob, patch, blockno, patchno);

    // Write Data Arrays
    int offset = WriteQDataArray(file, glob, patch);
    WriteIntDataArray(file, offset, "mpirank", glob->mpirank);
    offset += sizeof(uint32_t) + int_data_length;
    WriteIntDataArray(file, offset, "blockno", blockno);
    offset += sizeof(uint32_t) + int_data_length;
    WriteIntDataArray(file, offset, "patchno", patchno);
    
    // middle stuff
    fprintf(file, "\t\t\t</CellData>\n");
    fprintf(file, "\t\t\t<PointData>\n");
    fprintf(file, "\t\t\t</PointData>\n");
    fprintf(file, "\t\t</Piece>\n");
    fprintf(file, "\t</ImageData>\n");
    fprintf(file, "\t<AppendedData encoding=\"raw\">\n");
    fprintf(file, "\t\t_");

    //append data
    AppendQData(file, glob, patch);
    AppendIntData(file, glob->mpirank, int_data_length);
    AppendIntData(file, blockno, int_data_length);
    AppendIntData(file, patchno, int_data_length);


    //end file
    fprintf(file, "\n");
    fprintf(file, "\t</AppendedData>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
}

static void 
WriteVTPDFile(const char* basename, const fclaw_global_t *glob)
{
    char full_file_name[1024];
    snprintf(full_file_name, sizeof(full_file_name), "%s.vtpd", basename);

    FILE *file = fopen(full_file_name, "w");
    if (file == NULL) {
        perror("Could not open file");
        exit(1);
    }

    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<VTKFile type=\"vtkPartitionedDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    fprintf(file, "\t<vtkPartitionedDataSet>\n");

    for (int i = 0; i < glob->domain->global_num_patches; ++i) {
        char data_set_file[1024];
        snprintf(data_set_file, sizeof(data_set_file), "%s/%s_%d.vti", basename, basename, i);
        fprintf(file, "\t\t<DataSet index=\"%d\" file=\"%s\"/>\n", i, data_set_file);
    }

    fprintf(file, "\t</vtkPartitionedDataSet>\n");
    fprintf(file, "</VTKFile>\n");

    fclose(file);
}

void fclaw_clawpatch_output_vtpd_to_file (fclaw_global_t * glob, const char* basename)
{
    int rank = glob->mpirank;
    if (rank == 0) {
        WriteVTPDFile(basename, glob);
        fclaw_mkdir(basename);
    }
    sc_MPI_Barrier(glob->mpicomm);
    fclaw_global_iterate_patches(glob, WriteVTIFile,(void*) basename);
}
void fclaw_clawpatch_output_vtpd (fclaw_global_t * glob, int iframe)
{
    const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    char basename[BUFSIZ];
    snprintf (basename, BUFSIZ, "%s_frame_%04d", fclaw_opt->prefix, iframe);
    fclaw_clawpatch_output_vtpd_to_file(glob, basename);
}
