#ifndef FC2D_CUDA_PROFILER_H
#define FC2D_CUDA_PROFILER_H

#include <stdint.h>

class CudaTracer {
    public:
        CudaTracer(const char* name, int cid = 0);
        ~CudaTracer();
};

#ifdef FCLAW_ENABLE_DEBUG

#define PROFILE_CUDA(fname) CudaTracer uniq_name_using_macros__(fname);
#define PROFILE_CUDA_GROUP(fname, groupid) CudaTracer uniq_name_using_macros__(fname, groupid);

#else

#define PROFILE_CUDA(fname) ;
#define PROFILE_CUDA_GROUP(fname, groupid) ;

#endif

#endif
