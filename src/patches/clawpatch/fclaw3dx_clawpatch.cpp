/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#include <fclaw3dx_clawpatch.h>
#include <fclaw3dx_clawpatch.hpp>

#include <fclaw3dx_clawpatch_diagnostics.h>
#include <fclaw3dx_clawpatch_options.h>
#include <fclaw3dx_clawpatch_output_ascii.h> 
#include <fclaw3dx_clawpatch_output_vtk.h>
#include <fclaw3dx_clawpatch_fort.h>
#include <fclaw3dx_clawpatch_conservation.h>
#include <fclaw3dx_clawpatch_conservation_fort.h>
#include <fclaw3dx_clawpatch_transform.h>
#include <fclaw3dx_clawpatch_pillow.h>  

#include <fclaw3dx_clawpatch46_fort.h>

#include <_fclaw2d_to_fclaw3dx.h>

#include <fclaw2d_clawpatch.cpp>