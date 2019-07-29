#include <stdio.h>
#include <sys/time.h>

const int N_def (1 << 10);
const int threadsPerBlock = 32;
//const int blocksPerGrid = (N_def+threadsPerBlock-1) / threadsPerBlock;
const int blocksPerGrid = 10;

#define MAX(a,b) (a) > (b) ? (a) : (b)

__global__ void cuda_max(int N, double *a, double *c) 
{
    // __shared__ double localDot[threadsPerBlock];  /* Statically defined */
    extern __shared__ double localReduce[];
    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int localIndex = threadIdx.x;

    double localMax = a[ix];
    while (ix < N) 
    {
        localMax = MAX(localMax,a[ix]);
        ix += blockDim.x * gridDim.x;
    }
    
    /* Store local max computed by this thread */
    localReduce[localIndex] = localMax;
    
    /* Wait for all threads to get to this point */
    __syncthreads();

    /* Every block should compute max on the block */
    int i = blockDim.x/2;
    while (i != 0) 
    {
        if (localIndex < i)
        {
            localReduce[localIndex] = 
                  MAX(localReduce[localIndex],localReduce[localIndex+i]);
        }
        __syncthreads();
        i /= 2;
    }

    /* Each block stores local max */
    if (localIndex == 0)
        c[blockIdx.x] = localReduce[0];
}

double max_gpu(int N, double *a, 
               double *dev_a, 
               double *dev_partial_c)
{
    double   max, *partial_c;

    partial_c = (double*) malloc( blocksPerGrid*sizeof(double) );

    cudaMemcpy(dev_a, a, N*sizeof(double),
                              cudaMemcpyHostToDevice );

    dim3 block(threadsPerBlock);  /* Values defined in macros */
    dim3 grid(blocksPerGrid);     /* defined in macros, above */
    cuda_max<<<grid,block,threadsPerBlock*sizeof(double)>>>(N, dev_a,  
                                                            dev_partial_c );
    cudaDeviceSynchronize();
    cudaPeekAtLastError();


    /* copy the array 'c' back from the GPU to the CPU */
    cudaMemcpy( partial_c, dev_partial_c,
                      blocksPerGrid*sizeof(double),
                      cudaMemcpyDeviceToHost );

    /* Sum of block sums */
    max = partial_c[0];
    for (int i = 0; i < blocksPerGrid; i++) 
    {
        max = MAX(max,partial_c[i]);
    }

    free(partial_c);

    return max;
}


int main( void ) 
{
    double *a;
    double *dev_a, *dev_partial_c;
    double c_gpu;
    int N;
    double true_max = 3.14159;

    N = N_def;

    a = (double*) malloc( N*sizeof(double) );


    /* allocate the memory on the GPU */
    cudaMalloc((void**) &dev_a, N*sizeof(double));
    cudaMalloc((void**) &dev_partial_c, blocksPerGrid*sizeof(double) );


    /* Define vector a */
    for (int i = 0; i < N; i++) 
    {
        a[i] = 1.0;
    }
    a[N/3] = true_max;  /* Set the max here */

    c_gpu = max_gpu(N,a,dev_a,dev_partial_c);

    printf("%20s %10f\n","Maximum (GPU)", c_gpu);
    printf("%20s %10f\n","True maximum", true_max);

    cudaFree(dev_a);
    cudaFree(dev_partial_c);

    free(a);
}
