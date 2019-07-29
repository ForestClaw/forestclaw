#include <stdio.h>
#include <sys/time.h>

const int N_def (1 << 20);
const int threadsPerBlock = 32;
//const int blocksPerGrid = (N_def+threadsPerBlock-1) / threadsPerBlock;
const int blocksPerGrid = 1;


__global__ void cuda_dot(int N, double *a, double *b, double *c) 
{
    // __shared__ double localDot[threadsPerBlock];  /* Statically defined */
    extern __shared__ double localDot[];
    int ix = threadIdx.x + blockIdx.x * blockDim.x;
    int localIndex = threadIdx.x;

    double localSum = 0;
    while (ix < N) 
    {
        localSum += a[ix] * b[ix];  /* Reduction is here */
        ix += blockDim.x * gridDim.x;
    }
    
    /* Store sum computed by this thread */
    localDot[localIndex] = localSum;
    
    /* Wait for all threads to get to this point */
    __syncthreads();

    /* Every block should add up sum computed on  
       threads in the block */
    int i = blockDim.x/2;
    while (i != 0) 
    {
        if (localIndex < i)
        {
            localDot[localIndex] += localDot[localIndex + i];
        }
        __syncthreads();
        i /= 2;
    }

    /* Each block stores local dot product */
    if (localIndex == 0)
        c[blockIdx.x] = localDot[0];
}

double dot_gpu(int N, double *a, double *b,
               double *dev_a, double *dev_b, 
               double *dev_partial_c)
{
    double   dot, *partial_c;

    partial_c = (double*) malloc( blocksPerGrid*sizeof(double) );

    /* copy the arrays 'a' and 'b' to the GPU */
    cudaMemcpy(dev_a, a, N*sizeof(double),
                              cudaMemcpyHostToDevice );
    cudaMemcpy(dev_b, b, N*sizeof(double),
                              cudaMemcpyHostToDevice ); 

    dim3 block(threadsPerBlock);  /* Values defined in macros */
    dim3 grid(blocksPerGrid);     /* defined in macros, above */
    cuda_dot<<<grid,block,threadsPerBlock*sizeof(double)>>>(N, dev_a, dev_b, 
                                                            dev_partial_c );
    cudaDeviceSynchronize();
    cudaPeekAtLastError();


    /* copy the array 'c' back from the GPU to the CPU */
    cudaMemcpy( partial_c, dev_partial_c,
                      blocksPerGrid*sizeof(double),
                      cudaMemcpyDeviceToHost );

    /* Sum of block sums */
    dot = 0;
    for (int i = 0; i < blocksPerGrid; i++) 
    {
        dot += partial_c[i];
    }

    free(partial_c);

    return dot;
}

double dot_cpu(int n, double *a, double *b)
{
    double sum = 0;
    int i;

    for (i = 0; i < n; i++)
    {
        sum += a[i]*b[i];
    }
    return sum;
}

/* Compute a dot product */
int main( void ) 
{
    double *a, *b;
    double *dev_a, *dev_b, *dev_partial_c;
    double c_gpu;
    int N;

    N = N_def;

    a = (double*) malloc( N*sizeof(double) );
    b = (double*) malloc( N*sizeof(double) );


    /* allocate the memory on the GPU */
    cudaMalloc((void**) &dev_a, N*sizeof(double));
    cudaMalloc((void**) &dev_b, N*sizeof(double));
    cudaMalloc((void**) &dev_partial_c, blocksPerGrid*sizeof(double) );


    /* Define vectors a and b */
    for (int i = 0; i < N; i++) 
    {
        a[i] = 1.0;
        b[i] = 1.0;
    }

    /* GPU */
    c_gpu = dot_gpu(N,a,b,dev_a,dev_b,dev_partial_c);

    double s = N;   /* Sum of 1s */
    printf("%20s %10f\n","Dot product (GPU)", c_gpu);
    printf("%20s %10f\n","True dot product", s);

    /* free memory on the gpu side */
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_partial_c);

    free(a);
    free(b);
}
