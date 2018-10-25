/*
  Copyright (c) 2018 Carsten Burstedde, Donna Calhoun, Melody Shih, Scott Aiton, 
  Xinsheng Qin.
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


#include "../fc2d_cudaclaw_cuda.h"
#include "../fc2d_cudaclaw_check.cu"

#include <fclaw2d_global.h>
#include <fclaw_mpi.h>

void fc2d_cudaclaw_initialize_GPUs(fclaw2d_global_t *glob)
{
    cudaDeviceProp  prop;

    int mpirank, count, device_num;

    char name[MPI_MAX_PROCESSOR_NAME];
    int len;

    fclaw_global_essentialf("Block-size (FC2D_CUDACLAW_BLOCK_SIZE) set to %d\n",
                            FC2D_CUDACLAW_BLOCK_SIZE);            

    mpirank = glob->mpirank;

    CHECK(cudaGetDeviceCount(&count));

    device_num = mpirank % count;  

    CHECK(cudaSetDevice(device_num));


    /* Print out info */
    MPI_Get_processor_name(name, &len);

    fclaw_mpi_serialization_enter (glob);
    cudaGetDeviceProperties(&prop, device_num);
    printf("[fclaw] Rank %2d (%s) assigned to GPU %d  (%s)\n",mpirank, name, 
           device_num,prop.name); 
    fclaw_mpi_serialization_leave (glob);

    fflush(stdout);

    cudaDeviceReset();
}




