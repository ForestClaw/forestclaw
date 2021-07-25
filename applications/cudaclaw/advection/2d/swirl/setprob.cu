#include "../swirl_user.h"

#include <fc2d_cudaclaw_check.h>

void setprob()
{
    double period;
    FILE *f = fopen("setprob.data","r");
    fscanf(f,"%lf",&period);
    fclose(f);    

    CHECK(cudaMemcpyToSymbol(s_tperiod, &period, sizeof(double)));
}
