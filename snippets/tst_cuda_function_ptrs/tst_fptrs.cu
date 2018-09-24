#include <stdio.h>

typedef void (*fc2d_cuda_t)(float x,float *y);

__device__ void f1(float x, float *y)
{
    *y = 2*x;
    return;
}

__device__ void f2(float x, float *y)
{
    *y = -x;
    return;
}

__global__ void kernel(fc2d_cuda_t f,float x, float *y)
{
    f(x,y);
    return;
}

__device__ fc2d_cuda_t fptr = f1;

int main()
{
    fc2d_cuda_t h_f;

    cudaMemcpyFromSymbol(&h_f, fptr, sizeof(fc2d_cuda_t));

    float x = 5;
    float y;
    kernel<<<1,1>>>(h_f,x,&y);
    printf("x = %f; y = %f\n",x,y);

    return 0;
}

