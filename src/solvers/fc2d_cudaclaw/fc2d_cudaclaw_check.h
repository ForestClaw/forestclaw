#include <stdio.h>

#define FC2D_CUDACLAW_CHECK_H
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"ERROR : %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


