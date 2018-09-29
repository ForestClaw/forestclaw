#ifndef TST_FPTRS_H
#define TST_FPTRS_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef float (*fc2d_cuda_t)(float x);
typedef void (*fc2d_assign_cuda_ptr_t)(fc2d_cuda_t* h_f);

typedef struct fc2d_cuda_vt
{
    fc2d_cuda_t h_f;
} fc2d_cuda_vt_t;


#ifdef __cplusplus
}
#endif

#endif
