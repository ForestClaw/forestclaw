#include <fclaw_clawpatch_output_vtpd.h>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <vector>
#include <fclaw_global.h>
#include <fclaw_options.h>
#include <fclaw_clawpatch_options.h>
#include <fclaw_clawpatch.h>
using namespace std;
namespace
{
/**
 * @brief Write the dataset metadata information
 *
 * @param comm the communicator
 * @param file the file
 * @param domain the domain
 * @return int the offset in the file after the metadata information
 */
void 
WriteVTIFile(fclaw_domain_t * domain, fclaw_patch_t * patch,
     int blockno, int patchno, void *user)
{
	fclaw_global_iterate_t* g = (fclaw_global_iterate_t*) user;
	fclaw_global_t* glob = g->glob;
	char *base_name = (char*) g->user;

	const int globalno = glob->domain->blocks[blockno].num_patches_before + patchno;

	char file_name[BUFSIZ];
	snprintf(file_name, BUFSIZ, "%s/%s_%d.vti", base_name, base_name, globalno);

	ofstream file(file_name);

	int mx, my, mz, mbc;
	double xlower, ylower, zlower, dx, dy, dz;
	fclaw3d_clawpatch_grid_data(glob, patch, 
		&mx, &my, &mz, &mbc, 
		&xlower, &ylower, &zlower, 
		&dx, &dy, &dz);

	int start = (mx+2*mbc)*mbc + mbc;
	int stride_j = mx + 2 * mbc;
	int stride_k = stride_j * (my + 2 * mbc);

	int appended_data_stride = sizeof(uint32_t) + mx * my * mz * sizeof(double);

	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << endl;
	file << "\t<ImageData WholeExtent=\"0 " << mx << " 0 " << my << " 0 " << mz << "\" Origin=\"" << xlower << " " << ylower << " " << zlower << "\" Spacing=\"" << dx << " " << dy << " " << dz << "\" Direction=\"1 0 0 0 1 0 0 0 1\">" << endl;
	file << "\t\t<Piece Extent=\"0 " << mx << " 0 " << mx << " 0 " << mz << "\">" << endl;
	file << "\t\t\t<PointData>" << endl;
	file << "\t\t\t</PointData>" << endl;
	file << "\t\t\t<CellData>" << endl;
	int index = 0;
	for(int i=0; i<1; i++)
	{
		string name   = "meqn";
		double* q = fclaw_clawpatch_get_q(glob,patch);

		double min_val = numeric_limits<double>::infinity();
		double max_val = -numeric_limits<double>::infinity();

		for(int k = 0; k < mz; k++)
		for(int j = 0; j < my; j++)
		for(int i = 0; i < my; i++)
		{
			double val = q[i+mbc + (j+mbc) * stride_j + (k+mbc) * stride_k];
			min_val = min(min_val, val);
			max_val = max(max_val, val);
		}

		int offset = appended_data_stride * index;
		file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"appended\" RangeMin=\"" << min_val << "\" RangeMax=\"" << max_val << "\" offset=\"" << offset << "\">" << endl;
		file << "\t\t\t\t</DataArray>" << endl;
		index++;
	}
	file << "\t\t\t</CellData>" << endl;
	file << "\t\t</Piece>" << endl;
	file << "\t</ImageData>" << endl;
	file << "\t<AppendedData encoding=\"raw\">" << endl;
	file << "\t\t_";

	uint32_t size = mx * my * mz * sizeof(double);

	for(int i=0; i<1; i++)
	{
		string name   = "meqn";
		double* q = fclaw_clawpatch_get_q(glob,patch);

		file.write(reinterpret_cast<char *>(&size), sizeof(size));

		for(int k = 0; k < mz; k++)
		for(int j = 0; j < my; j++)
		for(int i = 0; i < my; i++)
		{
			double val = q[i+mbc + (j+mbc) * stride_j + (k+mbc) * stride_k];
			file.write(reinterpret_cast<const char *>(&val), sizeof(double));
		}
	}

	file << endl;
	file << "\t</AppendedData>" << endl;
	file << "</VTKFile>" << endl;

	file.close();
}
void WriteVTPDFile(const string &file_name, const fclaw_global_t *glob)
{
	ofstream file(file_name + ".vtpd");
	file << "<?xml version=\"1.0\"?>" << endl;
	file << "<VTKFile type=\"vtkPartitionedDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">" << endl;
	file << "\t<vtkPartitionedDataSet>" << endl;
	for (int i = 0; i < glob->domain->global_num_patches; i++) {
		file << "\t\t<DataSet index=\"" << i << "\" file=\"" << file_name << "/" << file_name << "_" << i << ".vti\"/>" << endl;
	}
	file << "\t</vtkPartitionedDataSet>" << endl;
	file << "</VTKFile>" << endl;
	file.close();
}
} // namespace
void fclaw_clawpatch_output_vtpd_to_file (fclaw_global_t * glob, const char* basename)
{
	int rank = glob->mpirank;
	if (rank == 0) {
		WriteVTPDFile(basename, glob);
		filesystem::create_directory(basename);
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
