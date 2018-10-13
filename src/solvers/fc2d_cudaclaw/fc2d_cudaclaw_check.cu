#include <stdio.h>

#define CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"ERROR : %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


