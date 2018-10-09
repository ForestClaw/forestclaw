/*
  Copyright (c) 2018 Carsten Burstedde, Donna Calhoun, Melody Shih, Scott Aiton, 
  Xinsheng Qin.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "../fc2d_cudaclaw_cuda.h" /* Only needed for external C code */

static
__device__ double minmod(double r)
{
    return fmax(0.0,fmin(1.0,r));
}

static
__device__ double superbee(double r)
{
    return fmax(0., fmax(fmin(1., 2.*r), fmin(2., r)));
}

static
__device__ double vanleer(double r)
{
    return (r + fabs(r)) / (1. + fabs(r));
}

static 
__device__ double monotinizedcentered(double r)
{
    double c = (1. + r)/2.;
    return fmax(0., fmin(c, fmin(2., 2.*r)));
}

static
__device__ double beamwarming(double r)
{
    return r;
}

__device__ double cudaclaw_limiter(int lim_choice, double r)
{
    switch(lim_choice)
    {
        case 1:
            return minmod(r);

        case 2:
            return superbee(r);

        case 3:
            return vanleer(r);

        case 4:
            return monotinizedcentered(r);

        case 5:
            return beamwarming(r);

        default:
            return 1;  /* No limiting */
    }
}
