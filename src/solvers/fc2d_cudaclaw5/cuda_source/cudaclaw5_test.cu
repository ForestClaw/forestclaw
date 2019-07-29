#include "../fc2d_cudaclaw5.h"

__global__ void donothing()
{
    /* Do nothing! */

    return;
}

__device__ int addem( int a, int b ) 
{
    return a + b;
}

__global__ void add( int a, int b, int *c ) 
{
    *c = a+b;
}


void cudaclaw5_test()
{
    int a,b,c;
    int *dev_c;

    donothing<<<1,1>>>();

    /* Allocate memory on the device */
    cudaMalloc( (void**)&dev_c, sizeof(int));

    a = 2;
    b = 7;
    c=2;
    add<<<1,1>>>(a, b, dev_c );   
	cudaError_t code = cudaPeekAtLastError();
    if(code!=cudaSuccess){
        printf("ERROR: %s\n",cudaGetErrorString(code));
    }


    cudaMemcpy( &c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);
 
    printf("Kernel result : %d + %d = %d\n",a,b,c);   

    cudaFree( dev_c);

    return;
}
