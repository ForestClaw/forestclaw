#include "tst_fptrs.h"
#include "user_f.h"

#include <stdio.h>

__global__ void kernel(fc2d_cuda_t f,float x, float *y)
{
    *y = f(x);
    return;
}

int main()
{
    float x;
    float y, *y_dev;
    fc2d_assign_cuda_ptr_t f_assign_user;

    fc2d_cuda_t h_f;

    /* User definitions */
    x = 5;
    f_assign_user = assign_cuda_ptr2;  /* Function pointer stored in a v-table */


    /* Code */
    f_assign_user(&h_f);    

    cudaMalloc((void**) &y_dev, sizeof(float));

    kernel<<<1,1>>>(h_f,x,y_dev);

    cudaMemcpy(&y, y_dev, sizeof(float), cudaMemcpyDeviceToHost);

    printf("x = %f; y = %f\n",x,y);

    return 0;
}

