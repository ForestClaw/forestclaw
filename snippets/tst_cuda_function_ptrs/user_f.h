#ifndef USER_F_H
#define USER_F_H

/* User defined functions */

#include "tst_fptrs.h"  /* Needed for type defs */

#ifdef __cplusplus
extern "C"
{
#endif

void assign_cuda_ptr1(fc2d_cuda_t* h_f);
void assign_cuda_ptr2(fc2d_cuda_t* h_f);

#ifdef __cplusplus
}
#endif

#endif


