#include "user_f.h"
#include <stdio.h>



__device__ float f1(float x)
{
    return 2*x;
}

__device__ fc2d_cuda_t fptr1 = f1;

void assign_cuda_ptr1(fc2d_cuda_t* h_f)
{
    cudaError_t ce = cudaMemcpyFromSymbol(h_f, fptr1, sizeof(fc2d_cuda_t));
    if(ce != cudaSuccess)
    {
        printf("ERROR: %s\n",cudaGetErrorString(ce));
        exit(0);
    }    
    else
    {
        printf("Success!\n");
    }
}

__device__ float f2(float x)
{
    return -x;
}

__device__ fc2d_cuda_t fptr2 = f2;

void assign_cuda_ptr2(fc2d_cuda_t* h_f)
{
    cudaError_t ce = cudaMemcpyFromSymbol(h_f, fptr2, sizeof(fc2d_cuda_t));
    if(ce != cudaSuccess)
    {
        printf("ERROR: %s\n",cudaGetErrorString(ce));
        exit(0);
    }      
    else
    {
        printf("Success!\n");
    }
  
}


