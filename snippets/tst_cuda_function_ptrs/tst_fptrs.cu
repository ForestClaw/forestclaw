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
    fc2d_cuda_vt_t vt;

    //fc2d_assign_cuda_ptr_t f_assign_user;

    fc2d_cuda_t h_f;

    /* User definitions (in swirl_user, for example) */
    x = 5;
    assign_cuda_ptr2(&vt.h_f);


    /* Code */
    cudaMalloc((void**) &y_dev, sizeof(float));

    kernel<<<1,1>>>(vt.h_f,x,y_dev);

    cudaMemcpy(&y, y_dev, sizeof(float), cudaMemcpyDeviceToHost);

    printf("x = %f; y = %f\n",x,y);

    return 0;
}

