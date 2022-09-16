#ifndef FC2D_CUDA_PROFILER_H
#define FC2D_CUDA_PROFILER_H

#include <stdint.h>

#ifdef FCLAW_ENABLE_DEBUG
#include <cuda_runtime_api.h>
#include <nvToolsExt.h>
#endif

class CudaTracer {
    public:
        CudaTracer(const char* name, int cid = 0);
        ~CudaTracer();
};

#ifdef FCLAW_ENABLE_DEBUG

#define CONCAT(a, b) CONCAT_INNER(a, b)
#define CONCAT_INNER(a, b) a ## b

#define UNIQUE_NAME(base) CONCAT(base, __COUNTER__)

#define PROFILE_CUDA(fname) CudaTracer UNIQUE_NAME(__tracer_)(fname);
#define PROFILE_CUDA_GROUP(fname, groupid) CudaTracer UNIQUE_NAME(__tracer_)(fname, groupid);

#else

#define PROFILE_CUDA(fname) ;
#define PROFILE_CUDA_GROUP(fname, groupid) ;

#endif

#endif
